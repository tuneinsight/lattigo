# DRLWE
The DRLWE package implements several ring-learning-with-errors-based Multiparty Homomorphic Encryption (MHE) primitives.
It provides generic interfaces for the local steps of the MHE-based Secure Multiparty Computation (MHE-MPC) protocol that are common between all the RLWE distributed schemes implemented in Lattigo (e.g., collective key generation).
The `dbfv` and `dckks` packages import the `drlwe` and provide scheme-specific functionalities (e.g., collective bootstrapping/refresh).

This package implements local operations only, hence does not assume nor provide any network-layer protocol implementation.
However, it provides a serialization methods for all relevant structures that implement the standard `encoding.BinaryMarshaller` and `encoding.BinaryUnmarshaller` interfaces (see [https://pkg.go.dev/encoding](https://pkg.go.dev/encoding)).

The MHE-MPC protocol implemented in Lattigo is based on the constructions described in ["Multiparty Homomorphic Encryption from Ring-Learning-with-Errors"](https://eprint.iacr.org/2020/304.pdf) by Mouchet et al. (2021), which is an RLWE instantiation of the MPC protocol described in ["Multiparty computation with low communication, computation and interaction via threshold FHE"](https://eprint.iacr.org/2011/613.pdf) by Asharov et al. (2012).

## MHE-MPC Protocol Overview

The protocol enables a group of N _input parties_ to compute a joint function over their private inputs under encryption and to provide a _receiver party_ with the result.
The protocol is generic and covers several system- and adversary-models:

**Peer-to-peer vs Cloud-assisted models**. The parties can execute the protocol in a peer-to-peer way or receive assistance from third-party server (which is also considered an adversary).

**Internal vs External receivers**. _Receiver parties_ are the intended recipients of the computation result and can be either _internal_ or _external_ depending or whether they are or are not input parties themselves, respectively. This distinction is important in practice, because external receivers do not need to be online (and even to be known) for the setup phase.

**Anytrust vs Full-threshold Access-structure**. As for many MPC protocols, the assumption on the worst-case number of corrupted parties can be reflected in the cryptographic access-control mechanism (the _access structure_). The implemented MHE-MPC protocol is "anytrust" (N-out-of-N-threshold) by default but can be relaxed to any positive threshold t-out-of-N (see Threshold Secret-Key Generation).

**Passive vs Active Adversaries**. The implemented MHE-MPC protocol is secure against passive adversaries, and can in theory be extended to active security by requiring the parties to produce proofs that their shares are correctly computed for every round. Note that those proof are not implemented in Lattigo.

An execution of the MHE-based MPC protocol has two phase, the Setup and the Evaluation phases, each of which comprises a number of sub-protocols as depicted below (the details of each protocols are provided later).

1. Setup Phase
    1. Secret Keys Generation        
    2. _[if t < N]_ Threshold Secret-Key Generation
    3. Collective Public Encryption-Key Generation
    4. Collective Public Evaluation-Key Generation
        1. Relinearization-Key
        2. Other required Switching-Keys
2. Evaluation Phase
   1. Input (Encryption)
   2. Circuit Evaluation
   3. Output phase  (Decryption)
      1. Collective Key-Switching
      2. Local Decryption


## MHE-MPC Protocol Steps Description
This section provides a description for each sub-protocol composing the MHE-MPC protocol and provide pointers to the relevant Lattigo types and methods.
This description is a first draft and is meant to be improved in the future.
For concrete code examples, see the `example/dbfv` and `example/drlwe` folders.
For a more formal exposition, see ["Multiparty Homomorphic Encryption from Ring-Learning-with-Errors"](https://eprint.iacr.org/2020/304.pdf) and [An Efficient Threshold Access-Structure for RLWE-Based Multiparty Homomorphic Encryption](https://eprint.iacr.org/2022/780).

The system model is abstracted by considering that the parties have access to a common public authenticated channel.
In the cloud-assisted setting, this public channel can be the helper cloud-server.
In the peer-to-peer setting, it could be a public broadcast channel.
We also assume that parties can communicate over private authenticated channels.

Several protocols require the parties to have access to common uniformly random polynomials (CRP), which are sampled from a common random string (CRS).
This CRS is implemented as an interface type `drlwe.CRS` that can be read from the parties as a part of the protocols (see below).
The `drlwe.CRS` can be implemented by a `utils.KeyedPRNG` type for which all parties use the same key.

### 1. Setup
In this phase, the parties generate the various keys that are required by the Evaluation phase.
Similarly to LSSS-based MPC protocols such as SPDZ, the setup phase does not depend on the input and can be pre-computed.
However, unlike LSSS-based MPC, the setup produces public-key material that can be re-used for an unlimited number of evaluation phases.


#### 1.i Secret Keys Generation
The parties generate their individual secret-keys locally by using a `rlwe.KeyGenerator`; this provides them with a `rlwe.SecretKey` type.
See [rlwe/keygen.go](../rlwe/keygen.go) further information on key-generation.

The _ideal secret-key_ is implicitly defined as the sum of all secret-keys.
Hence, this secret-key enforces an _N-out-N_ access structure which requires all the parties to collaborate in a ciphertext decryption (hence, tolerates N-1 dishonest parties).

#### 1.ii _[if t < N]_ Threshold Secret-Key Generation
For settings where an _N-out-N_ access structure is too restrictive (e.g., from an availability point of view), an optional Threshold Secret-Key Generation Protocol can be performed to enable _t-out-of-N_ access-structures (hence tolerating t-1 dishonest parties).
The idea of this protocol is to apply Shamir Secret Sharing to the _ideal secret-key_ in such a way that any group of _t_ parties can reconstruct it.
This is achieved by a single-round protocol where each party applies Shamir Secret-Sharing to its own share of the _ideal secret-key_.

We assume that each party is associated with a distinct `drlwe.ShamirPublicPoint` that is known to the other parties.

This protocol is implemented by the `drlwe.Thresholdizer` type and its steps are as follows:
- Each party generates a `drlwe.ShamirPolynomial` by using the `Thresholdizer.GenShamirPolynomial` method, then generates a share of type `drlwe.ShamirSecretShare` for each of the other parties' `ShamirPublicPoint` by using the `Thresholdizer.GenShamirSecretShare`.
- Each party privately sends the respective `ShamirSecretShare` to each of the other parties. 
- Each party aggregates all the `ShamirSecretShare`s it received using the `Thresholdizer.AggregateShares` method.

Each party stores its aggregated `ShamirSecretShare` for later use.

#### 1.iii Public Key Generation
The parties execute the collective public encryption-key generation protocol to obtain an encryption-key for the _ideal secret-key_.

The protocol is implemented by the `drlwe.CKGProtocol` type and its steps are as follows:
- Each party samples a common random polynomial (`drlwe.CKGCRP`) from the CRS by using the `CKGProtocol.SampleCRP` method.
- _[if t < N]_ Each party uses the `drlwe.Combiner.GenAdditiveShare` to obtain a t-out-of-t sharing and use the result as their secret-key in the next step.
- Each party generates a share (`drlwe.CKGShare`) from the CRP and their secret-key, by using the `CKGProtocol.GenShare` method.
- Each party discloses its share over the public channel. The shares are aggregated with the `CKGProtocol.AggregateShare` method.
- Each party can derive the public encryption-key (`rlwe.PublicKey`) by using the `CKGProtocol.GenPublicKey` method.

After the execution of this protocol, the parties have access to the collective public encryption-key, hence can provide their inputs to computations. 

#### 1.iv Evaluation-Key Generation
In order to evaluate circuits on the collectively-encrypted inputs, the parties must generate the switching-keys that correspond to the operations they wish to support.
The generation of a relinearization-key, that enables compact homomorphic multiplication, is described below (see `drlwe.RKGProtocol` below).
Additionally, and given that the circuit requires it, the parties can generate switching-keys to support rotations and other kinds of automorphisms (see `drlwe.RTGProtocol` below).

##### 1.iv.a Relinearization Key
This protocol provides the parties with a public relinearization-key (`rlwe.RelinearizationKey`) for the _ideal secret-key_. This public-key enables compact multiplications in RLWE schemes. This protocol is the only one that has two rounds.

The protocol is implemented by the  `drlwe.RKGProtocol` type and its steps are as follows:
- Each party samples a common random polynomial matrix (`drlwe.RKGCRP`) from the CRS by using the `RKGProtocol.SampleCRP` method.
- _[if t < N]_ Each party uses the `drlwe.Combiner.GenAdditiveShare` to obtain a t-out-of-t sharing and use the result as their secret-key in the next steps.
- Each party generates a share (`drlwe.RGKShare`) for the first protocol round by using the `RKGProtocol.GenShareRoundOne` method. This method also provides the party with an ephemeral secret-key (`rlwe.SecretKey`) that is required for the second round.
- Each party discloses its share for the first round over the public channel. The shares are aggregated with the `RKGProtocol.AggregateShare` method.
- Each party generates a share (also a `drlwe.RGKShare`) for the second protocol round by using the `RKGProtocol.GenShareRoundTwo` method.
- Each party discloses its share for the second round over the public channel. The shares are aggregated with the `RKGProtocol.AggregateShare` method.
- Each party can derive the public relinearization-key (`rlwe.RelinearizationKey`) by using the `RKGProtocol.GenRelinearizationKey` method.

#### 1.iv.b Rotation-keys and other Automorphisms
This protocol provides the parties with a public rotation-key (stored as `rlwe.SwitchingKey` types) for the _ideal secret-key_. One rotation-key enables one specific rotation on the ciphertexts' slots. The protocol can be repeated to generate the keys for multiple rotations.

The protocol is implemented by the  `drlwe.RTGProtocol` type and its steps are as follows:
- Each party samples a common random polynomial matrix (`drlwe.RTGCRP`) from the CRS by using the `RTGProtocol.SampleCRP` method.
- _[if t < N]_ Each party uses the `drlwe.Combiner.GenAdditiveShare` to obtain a t-out-of-t sharing and use the result as their secret-key in the next step.
- Each party generates a share (`drlwe.RTGShare`) by using `RTGProtocol.GenShare`. 
- Each party discloses its `drlwe.RTGShare` over the public channel. The shares are aggregated with the `RTGProtocol.AggregateShare` method.
- Each party can derive the public rotation-key (`rlwe.SwitchingKey`) from the final `RTGShare` by using the `RTGProtocol.AggregateShare` method.

### 2 Evaluation Phase 

#### 2.i Input step (Encryption Phase)

The parties provide their inputs for the computation during the Input Phase.
They use the collective encryption-key generated during the Setup Phase to encrypt their inputs, and disclose them on the public channel.
Since the collective encryption-key is a valid RLWE public encryption-key, it can be used directly with the single-party scheme.
Hence, the parties can use the `Encoder` and `Encryptor` interfaces of the desired encryption scheme (see [bfv.Encoder](../bfv/encoder.go), [bfv.Encryptor](../bfv/encryptor.go), [ckks.Encoder](../ckks/encoder.go) and [ckks.Encryptor](../ckks/encryptor.go)).


#### 2.ii Circuit Evaluation step
The computation of the desired function is performed homomorphically during the Evaluation Phase.
The step can be performed by the parties themselves or can be outsourced to a cloud-server. 
Since the ciphertexts in the multiparty schemes are valid ciphertexts for the single-party ones, the homomorphic operation of the latter can be used directly (see [bfv.Evaluator](../bfv/evaluator.go) and [ckks.Evaluator](../ckks/evaluator.go)).


#### 2.iii Output step
The receiver(s) obtain their outputs through the final Output Phase, which aim is to decrypt the ciphertexts resulting from the Evaluation Phase.
It is a two-steps process with an optional pre-processing step in case of the t-out-of-N access-structure.
In the first step, Collective Key-Switching the parties re-encrypt the desired ciphertext under the receiver's secret-key.
The second step is the local decryption of this re-encrypted ciphertext by the receiver.

#### 2.iii.a Collective Key-Switching
The parties perform a re-encryption of the desired ciphertext(s) from being encrypted under the _ideal secret-key_ to being encrypted under the receiver's secret-key.
There are two instantiations of the Collective Key-Switching protocol:
- Collective Key-Switching (CKS), implemented as the `drlwe.CKSProtocol` interface, enables the parties to switch from their _ideal secret-key_ _s_ to another _ideal secret-key_ _s'_ when s' is collectively known to the parties. In the case where _s' = 0_, this is equivalent to a decryption protocol which can be used when the receiver is one of the input-parties. 
- Collective Public-Key Switching (PCKS), implemented as the `drlwe.PCKSProtocol` interface, enables parties to switch from their _ideal secret-key_ _s_ to an arbitrary key _s'_ when provided with a public encryption-key for _s'_. Hence, this enables key-switching to a secret-key that is not known to the input parties, which further enables external receivers.

While both protocol variants have slightly different local operation, their steps are the same and as follows:
- Each party generates a share (of type `drlwe.CKSShare` or `drlwe.PCKSShare`) with the `drlwe.(P)CKSProtocol.GenShare` method. This requires its own secret-key (a `rlwe.SecretKey`) as well as the destination key: its own share of the destination key (a `rlwe.SecretKey`) in CKS or the destination public-key (a `rlwe.PublicKey`) in PCKS.
- Each party discloses its `drlwe.CKSShare` over the public channel. The shares are aggregated with the `(P)CKSProtocol.AggregateShares` method.
- From the aggregated `drlwe.CKSShare`, any party can derive the ciphertext re-encrypted under _s'_ by using the `(P)CKSProtocol.KeySwitch` method.

#### 2.iii.c Decryption
Once the receivers have obtained the ciphertext re-encrypted under their respective keys, they can use the usual decryption algorithm of the single-party scheme to obtain the plaintext result (see [bfv.Decryptor](../bfv/decryptor.go) and [ckks.Decryptor](../ckks/decryptor.go)).

