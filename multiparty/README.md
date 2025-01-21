# Multiparty Schemes in Lattigo

The `multiparty` package implements several Multiparty Homomorphic Encryption (MHE) primitives based on Ring-Learning-with-Errors (RLWE).
It provides the implementation of two core schemes:

1. A $N$-out-of-$N$-threshold scheme
2. A $t$-out-of-$N$-threshold scheme

We provide more informations about these two core schemes below. Moreover, The `multiparty/mpbgv` and `multiparty/mpckks` packages provide scheme-specific functionalities (e.g., interactive bootstrapping) by implementing **threshold** versions of the single-party BFV/BGV and CKKS cryptosystems
found in the `schemes` package.
Note that, as for the single party schemes, most of the operations are generic and can be handled by the core schemes in the `multiparty` package.

The `multiparty` package and its sub-packages provide the **local steps** of the MHE-based Secure Multiparty Computation protocol (MHE-MPC), as described in ["Multiparty Homomorphic Encryption from Ring-Learning-with-Errors"](https://eprint.iacr.org/2020/304.pdf) by Mouchet et al. (2021) (which is an RLWE instantiation of the MPC protocol described in ["Multiparty computation with low communication, computation and interaction via threshold FHE"](https://eprint.iacr.org/2011/613.pdf) by Asharov et al. (2012)).

Note that this package implements local operations only, hence does not assume or provide any network-layer protocol implementation.
However:
- All relevant types provide serialization methods implementing the `encoding.BinaryMarshaller` and `encoding.BinaryUnmarshaller` interfaces (see [https://pkg.go.dev/encoding](https://pkg.go.dev/encoding)) as well as the `io.WriterTo` and `io.ReaderFrom` interfaces (see [https://pkg.go.dev/io](https://pkg.go.dev/io)).
- The last section of this README provides a detailed overview of the MHE-MPC protocol, its different instantiations, and maps the protocol steps to the relevant Lattigo types and methods.
- The `examples/multiparty` folder contains example applications for simple MPC tasks. These examples are running all the parties in the same process, but demonstrate the use of the multiparty schemes in the MHE-MPC protocol.

## The $N$-out-of-$N$-Threshold Scheme

Conceptually, the $N$-out-of-$N$-threshold scheme exploits the linearity of RLWE encryption to distribute the secret-key among $N$ parties. More specifically, the core cryptographic operation of (single-party) RLWE-based scheme is to compute functions of the form: $$F(a,s) = as+e$$
over a ring $R$ where $a \in R$ is public, $s \in R$ is the secret-key of the scheme and $e \in R$ is a small ring element (sampled fresh for each function). For example, notice that generating an RLWE public-key corresponds to exactly this operation. The $N$-out-of-$N$-threshold scheme consists in splitting the secret-key $s$ into $N$ additive shares such that $s=\sum^N_{i=1} s_i$ and that $s_i$ is held by party $i$. 
In this way, any secret-key operation (especially, decryption) requires the collaboration between **all** the $N$ parties.

Exploiting the linear structure of $s$ and the (almost) linearity of $F$, the latter can be computed piece-wise by the parties, as: 

$$h_i = F(a, s_i) = as_i + e_i.$$

Then $F(a,s)$ can be computed as $h=\sum^N_{i=1}h_i$. The output $h$ is the desired function, only with larger error. Moreover, the following features are relevant:

- Since $F(a, s_i)$ has the form of a RLWE sample, it can be publicly disclosed for aggregation. For example, the parties could use a helper for aggregating their shares $h_i$.
- Since the secret-key shares $s_i$ are random and never disclosed, they can be sampled locally by the parties (as normal RLWE secret-keys). Hence, no trusted dealer or private channels are necessary between the parties.
- Obtaining protocols for threshold generation of evaluation-keys and threshold decryption only requires to slightly adapt the function $F$ above. The overall linearity argument remains the same.
- The threshold decryption protocol can be generalized into a re-encryption protocol (where the actual decryption corresponds to re-encrypting towards the $s=0$ key). It is possible to re-encrypt towards a known public-key (for receivers external to the set of parties).
- The threshold decryption- and key-switching protocols require *smudging* parameter implementing the noise-flooding technique. This is a statistical security parameter ensuring that no secret value leaks to the decryption receiver through the error. See the SECURITY.md for more discussion on this aspect.

### Implementation

For each secret-key operation in the single party RLWE scheme, the  `multiparty` package provides a "protocol" type implementing the local operations of the parties in the threshold scheme. More specifically, each protocol type provides a `GenShare(*rlwe.SecretKey, [...])` method for generating a party's share (i.e., computing $F(a, s_i)$ above) and a `AggregateShares(share1, share2 [...])` method to aggregate the shares. Moreover:

- For some protocols, $a$ is a *common random polynomial* (CRP) sampled from a common random string (CRS). In this case, the protocol types provide a `SampleCRP(crs CRS)` method to obtain $a$.
- The protocol types also provide method for converting the final aggregated share into the protocol's output. E.g., the `PublicKeyGenProtocol` provides a `GenPublicKey(aggshare PublicKeyGenShare, [...])` method.

## The $t$-out-of-$N$-Threshold Scheme

There might be settings where an $N$-out-of-$N$-threshold access-structure is too restrictive. For example, when $N$ is large, the probability of a single party being down at a given time increases. In cases where it can be assumed that the adversary cannot corrupt more than $t-1$ out of the $N$ parties, the $t$-out-of-$N$-threshold scheme can be employed to provide better liveness guarantees. More specifically, this scheme ensures that secret-key operations can be performed by any group of at least $t$ parties.

Lattigo provides an implementation of the RLWE-based $t$-out-of-$N$-threshold scheme described in Mouchet et al.'s paper [An Efficient Threshold Access-Structure for RLWE-Based Multiparty Homomorphic Encryption](https://eprint.iacr.org/2022/780). Similarly to many threshold schemes, it relies on Shamir Secret Sharing to distribute the secret-key of the scheme. This is, the secret-key of the scheme is now of the form:

$$S(X) = s + s_1 X + s_2 X^2 + ... + s_t X^{t-1},$$

i.e., a degree-$(t-1)$ polynomial in R[X] for which $s = S(0)$, and party $i$'s secret-key shares is distributed as $S(\alpha_i)$ for $(\alpha_1, \alpha_2, ... \alpha_N)$ $N$ distinct elements of $R$ forming an exceptional sequence.
Then, observe that $s$ can be reconstructed from any set of $t$ shares via Lagrange interpolation. For example, assuming reconstruction from the first $t$ shares:

$$s = \sum^t_{i=1} S(\alpha_i) \cdot \prod^t_{\substack{j=1\\ i \neq j}} \frac{\alpha_j}{\alpha_j - \alpha_i} = \sum^t_{i=1} S(\alpha_i) \cdot l_i$$

Hence, the structure of $s$ is still linear, and we can again compute the $F$ function above piece-wise. However, the fact that $F(a,s)$ is only **approximately** linear in $s$ poses some challenge in performing the Shamir reconstruction in the "usual" way. 
More specifically, we cannot compute $h_i = F(a, S(\alpha_i))$ locally and then combine the shares as $\sum^t_{i=1} h_i \cdot l_i$. This is because $l_i$ is a large $R$ element multiplying it with $S(0)a+e_i$ would result in a large error $e_i \cdot l_i$.

The scheme of Mouchet et al. circumvents this issue by directly evaluating $h_i=F(a, S(\alpha_i) \cdot l_i)$ locally. Then the combination of the share is back to being a simple summation over $t$ shares: $h =\sum^t_{i=1} h_i$. This simple trick enables a very efficient and usable $t$-out-of-$N$ scheme:

- $S$ can be generated non-interactively and without a trusted dealer by having each party generating a random degree-$(t-1)$ polynomial $S_i$ with $S_i(0) = s_i$, and by implicitly take $S=\sum^N_{i=1} S_i$. Observe then that $s = S(0) = \sum^N_{i=1} s_i$, which matches the $N$-out-of-$N$-threshold case.
- Then, party $i$ can obtain its share $S(\alpha_i)$ by:
    1. having each party $j$ send $S_j(\alpha_i)$ to party $i$ (via a **private** channel),
    2. having party $i$ compute $S(\alpha_i) = \sum^N_{j=1} S_j(\alpha_i)$.
- The above protocol is a single-round protocol, and state each party has to keep is then a single ring element $S(\alpha_i)$. 
- When instantiated as above, the $t$-out-of-$N$-threshold scheme consists in a direct **extension** of the $N$-out-of-$N$-threshold scheme where:
    1. The parties operate a *re-sharing* of their secret-key $N$-out-of-$N$ secret-key share using the $t$-out-of-$N$ Shamir Secret Sharing scheme. 
    2. The parties perform the secret-key operations (i.e., the protocols) in the same way as in the $N$-out-of-$N$-threshold scheme, yet among $t$ parties only and with $S(\alpha_i)\cdot l_i$ instead of $s_i$.

However, the scheme has the downside of requiring to know set of parties participating to a given secret-key operation (i.e., evaluation of $F$). This is because evaluating $F(a, S(\alpha_i) \cdot l_i$ requires each party $i$ to compute the Lagrange coefficient $l_i$, which depends on the participating set. 
Another downside of this scheme is that it requires a round of private, pairwise message exchanges between the parties before the scheme can be used in the $t$-out-of-$N$ regime.

### Implementation

Thanks to the $t$-out-of-$N$-threshold scheme being a direct extension of the $N$-out-of-$N$-threshold scheme (see the discussion above), the implementation of the former consist of two new types: `Thresholdizer` and `Combiner`.

The `Thresholdizer` type implements the secret-key generation and re-sharing steps. This type corresponds to part 1. of the extension as described above. More specifically: 
- `GenShamirPolynomial(threshold int, sk *rlwe.SecretKey)` generates the random degree-$(t-1)$ polynomial with `sk` as constant coefficient (i.e., $S_i$ above).
- `GenShamirSecretShare(recipient ShamirPublicPoint, [...])` generates re-sharing for a given recipient by evaluating the secret polynomial (i.e., $S_i(\alpha_j)$ above).
- `AggregateShares(share1, share2 ShamirSecretShare, [...])` aggregates two received shares (i.e., one addition step in computing $S(\alpha_i)$ above).

The `Combiner` type lets parties obtain $t$-out-of-$t$ additive shares from their $t$-out-of-$N$ Shamir shares. This type corresponds to part 2. of the extension as described above, and is called as a pre-processing before any secret-key operation performed in the $t$-out-of-$N$ regime. More specifically, the `Combiner.GenAdditiveShare` takes as input the $t$-out-of-$N$-threshold secret-share of the party ($S(\alpha_i)$ above) along with the set $L=\{\alpha_1, ..., \alpha_t\}$ of the $t$ parties participating to the protocol, and computes:

$$S(\alpha_i) \cdot l_i = S(\alpha_i) \cdot \prod_{\substack{\alpha_j \in L\\ \alpha_j \neq \alpha_i}} \frac{\alpha_j}{\alpha_j - \alpha_i}.$$

Hence, from the share output by `GenAdditiveShare`, the usual protocols described for the $N$-out-of-$N$-threshold setting (see the previous section) can be used, yet with $N=t$.

## MHE-MPC Protocol Overview

The protocol enables a group of N _input parties_ to compute a joint function over their private inputs under encryption and to provide a _receiver party_ with the result.
The protocol is generic and covers several system- and adversary-models:

**Peer-to-peer vs Cloud-assisted models**. The parties can execute the protocol in a peer-to-peer way or receive assistance from a third-party server (which is also considered an adversary).

**Internal vs External receivers**. _Receiver parties_ are the intended recipients of the computation result and can be either _internal_ or _external_ depending or whether they are or not input parties themselves, respectively. This distinction is important in practice, because external receivers do not need to be online (and even to be known) during the setup phase.

**Anytrust vs Full-threshold Access-structure**. As for many MPC protocols, the assumption on the worst-case number of corrupted parties can be mapped in the cryptographic access-control mechanism (the _access structure_). The implemented MHE-MPC protocol is "anytrust" (N-out-of-N-threshold) by default, but can be relaxed to any positive threshold t-out-of-N (see Threshold Secret-Key Generation).

**Passive vs Active Adversaries**. The implemented MHE-MPC protocol is secure against passive adversaries, and can in theory be extended to active security by requiring the parties to produce proofs that their shares are correctly computed for every round. Note that those proofs are not implemented in Lattigo.

An execution of the MHE-based MPC protocol has two phases: the Setup phase and the Evaluation phase, each of which comprises a number of sub-protocols as depicted below (the details of each protocols are provided later).

1. Setup Phase
    1. Secret Keys Generation        
    2. _[if t < N]_ Threshold Secret-Key Generation
    3. Collective Public Encryption-Key Generation
    4. Collective Public Evaluation-Key Generation
        1. Relinearization-Key
        2. Galois-Keys
        3. Generic Evaluation-Keys
2. Evaluation Phase
   1. Input (Encryption)
   2. Circuit Evaluation
   3. Output phase (Decryption)
      1. Collective Key-Switching
      2. Local Decryption


## MHE-MPC Protocol Steps Description

This section provides a description for each sub-protocol of the MHE-MPC protocol and provides pointers to the relevant Lattigo types and methods.
This description is a first draft and will evolve in the future.
For concrete code examples, see the `example/multiparty` folders.
For a more formal exposition, see ["Multiparty Homomorphic Encryption from Ring-Learning-with-Errors"](https://eprint.iacr.org/2020/304.pdf) and [An Efficient Threshold Access-Structure for RLWE-Based Multiparty Homomorphic Encryption](https://eprint.iacr.org/2022/780).

The system model is abstracted by considering that the parties have access to a common public authenticated channel.
In the cloud-assisted setting, this public channel can be the helper cloud-server.
In the peer-to-peer setting, it could be a public broadcast channel.
We also assume that parties can communicate over private authenticated channels.

Several protocols require the parties to have access to common uniformly random polynomials (CRP), which are sampled from a common random string (CRS).
This CRS is implemented as an interface type `multiparty.CRS` that can be read by the parties as a part of the protocols (see below).
The `multiparty.CRS` can be implemented by a `utils.KeyedPRNG` type for which all parties use the same key.

### 1. Setup
In this phase, the parties generate the various keys that are required by the Evaluation phase.
Similarly to LSSS-based MPC protocols such as SPDZ, the setup phase does not depend on the input and can be pre-computed.
However, unlike LSSS-based MPC, the setup produces public-keys that can be re-used for an arbitrary number of evaluation phases.

#### 1.i Secret Keys Generation
The parties generate their individual secret-keys locally by using a `rlwe.KeyGenerator`; this provides them with a `rlwe.SecretKey` type.
See [core/rlwe/keygenerator.go](../core/rlwe/keygenerator.go) for further information on key-generation.

The _ideal secret-key_ is implicitly defined as the sum of all secret-keys.
Hence, this secret-key enforces an _N-out-N_ access structure which requires all the parties to collaborate in a ciphertext decryption and thus tolerates N-1 dishonest parties.

#### 1.ii _[if t < N]_ Threshold Secret-Key Generation
For settings where an _N-out-N_ access structure is too restrictive (e.g., from an availability point of view), an optional Threshold Secret-Key Generation Protocol can be run to enable _t-out-of-N_ access-structures (hence tolerating t-1 dishonest parties).
The idea of this protocol is to apply Shamir Secret Sharing to the _ideal secret-key_ in such a way that any group of _t_ parties can reconstruct it.
This is achieved by a single-round protocol where each party applies Shamir Secret-Sharing to its own share of the _ideal secret-key_.

We assume that each party is associated with a distinct `multiparty.ShamirPublicPoint` that is known to the other parties.

This protocol is implemented by the `multiparty.Thresholdizer` type and its steps are the following:
- Each party generates a `multiparty.ShamirPolynomial` by using the `Thresholdizer.GenShamirPolynomial` method, then generates a share of type `multiparty.ShamirSecretShare` for each of the other parties' `ShamirPublicPoint` by using the `Thresholdizer.GenShamirSecretShare`.
- Each party privately sends the respective `ShamirSecretShare` to each of the other parties. 
- Each party aggregates all the `ShamirSecretShare`s it received using the `Thresholdizer.AggregateShares` method.

Each party stores its aggregated `ShamirSecretShare` for later use.

#### 1.iii Public Key Generation
The parties execute the collective public encryption-key generation protocol to obtain an encryption-key for the _ideal secret-key_.

The protocol is implemented by the `multiparty.PublicKeyGenProtocol` type and its steps are as follows:
- Each party samples a common random polynomial (`multiparty.PublicKeyGenCRP`) from the CRS by using the `PublicKeyGenProtocol.SampleCRP` method.
- _[if t < N]_ Each party uses the `multiparty.Combiner.GenAdditiveShare` to obtain a t-out-of-t sharing and uses the result as its secret-key in the next step.
- Each party generates a share (`multiparty.PublicKeyGenShare`) from the CRP and their secret-key, by using the `PublicKeyGenProtocol.GenShare` method.
- Each party discloses its share over the public channel. The shares are aggregated with the `PublicKeyGenProtocol.AggregateShares` method.
- Each party can derive the public encryption-key (`rlwe.PublicKey`) by using the `PublicKeyGenProtocol.GenPublicKey` method.

After the execution of this protocol, the parties have access to the collective public encryption-key, hence can provide their inputs to computations. 

#### 1.iv Evaluation-Key Generation
In order to evaluate circuits on the collectively-encrypted inputs, the parties must generate the evaluation-keys that correspond to the operations they wish to support.
The generation of a relinearization-key, which enables compact homomorphic multiplication, is described below (see `multiparty.RelinearizationKeyGenProtocol`).
Additionally, and given that the circuit requires it, the parties can generate evaluation-keys to support rotations and other kinds of Galois automorphisms (see `multiparty.GaloisKeyGenProtocol` below).
Finally, it is possible to generate generic evaluation-keys to homomorphically re-encrypt a ciphertext from a secret-key to another (see `multiparty.EvaluationKeyGenProtocol`).

##### 1.iv.a Relinearization Key
This protocol provides the parties with a public relinearization-key (`rlwe.RelinearizationKey`) for the _ideal secret-key_. This public-key enables compact multiplications in RLWE schemes. Out of the described protocols in this package, this is the only two-round protocol.

The protocol is implemented by the  `multiparty.RelinearizationKeyGenProtocol` type and its steps are as follows:
- Each party samples a common random polynomial matrix (`multiparty.RelinearizationKeyGenCRP`) from the CRS by using the `RelinearizationKeyGenProtocol.SampleCRP` method.
- _[if t < N]_ Each party uses the `multiparty.Combiner.GenAdditiveShare` to obtain a t-out-of-t sharing and use the result as their secret-key in the next steps.
- Each party generates a share (`multiparty.RelinearizationKeyGenShare`) for the first protocol round by using the `RelinearizationKeyGenProtocol.GenShareRoundOne` method. This method also provides the party with an ephemeral secret-key (`rlwe.SecretKey`), which is required for the second round.
- Each party discloses its share for the first round over the public channel. The shares are aggregated with the `RelinearizationKeyGenProtocol.AggregateShares` method.
- Each party generates a share (also a `multiparty.RelinearizationKeyGenShare`) for the second protocol round by using the `RelinearizationKeyGenProtocol.GenShareRoundTwo` method.
- Each party discloses its share for the second round over the public channel. The shares are aggregated with the `RelinearizationKeyGenProtocol.AggregateShares` method.
- Each party can derive the public relinearization-key (`rlwe.RelinearizationKey`) by using the `RelinearizationKeyGenProtocol.GenRelinearizationKey` method.

##### 1.iv.b Galois Keys
This protocol provides the parties with a public Galois-key (stored as `rlwe.GaloisKey` types) for the _ideal secret-key_. One Galois-key enables one specific Galois automorphism on the ciphertexts' slots. The protocol can be repeated to generate the keys for multiple automorphisms.

The protocol is implemented by the `multiparty.GaloisKeyGenProtocol` type and its steps are as follows:
- Each party samples a common random polynomial matrix (`multiparty.GaloisKeyGenCRP`) from the CRS by using the `GaloisKeyGenProtocol.SampleCRP` method.
- _[if t < N]_ Each party uses the `multiparty.Combiner.GenAdditiveShare` to obtain a t-out-of-t sharing and uses the result as its secret-key in the next step.
- Each party generates a share (`multiparty.GaloisKeyGenShare`) by using `GaloisKeyGenProtocol.GenShare`. 
- Each party discloses its `multiparty.GaloisKeyGenShare` over the public channel. The shares are aggregated with the `GaloisKeyGenProtocol.AggregateShares` method.
- Each party can derive the public Galois-key (`rlwe.GaloisKey`) from the final `GaloisKeyGenShare` by using the `GaloisKeyGenProtocol.AggregateShares` method.

##### 1.iv.c Other Evaluation Keys
This protocol provides the parties with a generic public Evaluation-key (stored as `rlwe.EvaluationKey` types) for the _ideal secret-key_. One Evaluation-key enables one specific public re-encryption from one key to another.

The protocol is implemented by the  `multiparty.EvaluationKeyGenProtocol` type and its steps are as follows:
- Each party samples a common random polynomial matrix (`multiparty.EvaluationKeyGenCRP`) from the CRS by using the `EvaluationKeyGenProtocol.SampleCRP` method.
- _[if t < N]_ Each party uses the `multiparty.Combiner.GenAdditiveShare` to obtain a t-out-of-t sharing and uses the result as its secret-key in the next step.
- Each party generates a share (`multiparty.EvaluationKeyGenShare`) by using `EvaluationKeyGenProtocol.GenShare`. 
- Each party discloses its `multiparty.EvaluationKeyGenShare` over the public channel. The shares are aggregated with the `EvaluationKeyGenProtocol.AggregateShares` method.
- Each party can derive the public Evaluation-key (`rlwe.EvaluationKey`) from the final `EvaluationKeyGenShare` by using the `EvaluationKeyGenProtocol.AggregateShares` method.

### 2 Evaluation Phase 

#### 2.i Input step (Encryption Phase)

The parties provide their inputs for the computation during the Input Phase.
They use the collective encryption-key generated during the Setup Phase to encrypt their inputs, and send them through the public channel.
Since the collective encryption-key is a valid RLWE public encryption-key, it can be used directly with the single-party scheme.
Hence, the parties can use the `Encoder` and `Encryptor` interfaces of the desired encryption scheme (see [bgv.Encoder](../schemes/bgv/encoder.go), [ckks.Encoder](../schemes/ckks/encoder.go) and [rlwe.Encryptor](../core/rlwe/encryptor.go)).

#### 2.ii Circuit Evaluation step
The computation of the desired function is performed homomorphically during the Evaluation Phase.
The step can be performed by the parties themselves or can be outsourced to a cloud-server. 
Since the ciphertexts in the multiparty schemes are valid ciphertexts for the single-party ones, the homomorphic operation of the latter can be used directly (see [bgv.Evaluator](../schemes/bgv/evaluator.go) and [ckks.Evaluator](../schemes/ckks/evaluator.go)).

#### 2.iii Output step
The receiver(s) obtain their outputs through the final Output Phase, whose aim is to decrypt the ciphertexts resulting from the Evaluation Phase.
It is a two-step process with an optional pre-processing step when using the t-out-of-N access structure.
In the first step, Collective Key-Switching, the parties re-encrypt the desired ciphertext under the receiver's secret-key.
The second step is the local decryption of this re-encrypted ciphertext by the receiver.

##### 2.iii.a Collective Key-Switching
The parties perform a re-encryption of the desired ciphertext(s) from being encrypted under the _ideal secret-key_ to being encrypted under the receiver's secret-key.
There are two instantiations of the Collective Key-Switching protocol:
- Collective Key-Switching (KeySwitch), implemented as the `multiparty.KeySwitchProtocol` interface: it enables the parties to switch from their _ideal secret-key_ _s_ to another _ideal secret-key_ _s'_ when s' is collectively known by the parties. In the case where _s' = 0_, this is equivalent to a collective decryption protocol that can be used when the receiver is one of the input-parties. 
- Collective Public-Key Switching (PublicKeySwitch), implemented as the `multiparty.PublicKeySwitchProtocol` interface, enables parties to switch from their _ideal secret-key_ _s_ to an arbitrary key _s'_ when provided with a public encryption-key for _s'_. Hence, this enables key-switching to a secret-key that is not known to the input parties, which enables external receivers.

While both protocol variants have slightly different local operations, their steps are the same:
- Each party generates a share (of type `multiparty.KeySwitchShare` or `multiparty.PublicKeySwitchShare`) with the `multiparty.(Public)KeySwitchProtocol.GenShare` method. This requires its own secret-key (a `rlwe.SecretKey`) as well as the destination key: its own share of the destination key (a `rlwe.SecretKey`) in KeySwitch or the destination public-key (a `rlwe.PublicKey`) in PublicKeySwitch.
- Each party discloses its `multiparty.KeySwitchShare` over the public channel. The shares are aggregated with the `(Public)KeySwitchProtocol.AggregateShares` method.
- From the aggregated `multiparty.KeySwitchShare`, any party can derive the ciphertext re-encrypted under _s'_ by using the `(Public)KeySwitchProtocol.KeySwitch` method.

##### 2.iii.b Decryption
Once the receivers have obtained the ciphertext re-encrypted under their respective keys, they can use the usual decryption algorithm of the single-party scheme to obtain the plaintext result (see [rlwe.Decryptor](../core/rlwe/decryptor.go).

## References

1. Multiparty Homomorphic Encryption from Ring-Learning-With-Errors (<https://eprint.iacr.org/2020/304>)
2. An Efficient Threshold Access-Structure for RLWE-Based Multiparty Homomorphic Encryption (<https://eprint.iacr.org/2022/780>)
