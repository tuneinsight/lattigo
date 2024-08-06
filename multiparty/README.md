# Multi-Party

The `multiparty` package implements several Multiparty Homomorphic Encryption (MHE) primitives based on Ring-Learning-with-Errors (RLWE).
It provides generic interfaces for the local steps of the MHE-based Secure Multiparty Computation (MHE-MPC) protocol that are common across all the RLWE distributed schemes implemented in Lattigo (e.g., collective key generation).
The `multiparty/mpbgv` and `multiparty/mpckks` packages provide scheme-specific functionalities (e.g., interactive bootstrapping) by implementing **threshold** versions of the single-party BFV/BGV and CKKS cryptosystems
found in the `schemes` package.

This package implements local operations only, hence does not assume or provide any network-layer protocol implementation.
However, it provides serialization methods for all relevant structures that implement the standard `encoding.BinaryMarshaller` and `encoding.BinaryUnmarshaller` interfaces (see [https://pkg.go.dev/encoding](https://pkg.go.dev/encoding)) as well as the `io.WriterTo` and `io.ReaderFrom` interfaces (see [https://pkg.go.dev/encoding](https://pkg.go.dev/io)).

The MHE-MPC protocol implemented in Lattigo is based on the constructions described in ["Multiparty Homomorphic Encryption from Ring-Learning-with-Errors"](https://eprint.iacr.org/2020/304.pdf) by Mouchet et al. (2021), which is an RLWE instantiation of the MPC protocol described in ["Multiparty computation with low communication, computation and interaction via threshold FHE"](https://eprint.iacr.org/2011/613.pdf) by Asharov et al. (2012).

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
