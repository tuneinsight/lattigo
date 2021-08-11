# DRLWE
The DRLWE package is a collection of protocols to provide Multiparty Homomorphic Encryption (MHE) based on the ring-learning-with-errors.
It provides generic interfaces that implement the steps of a the MHE-based Secure Multiparty Computation (MPC) protocol that are common between all the RLWE distributed schemes implemented in Lattigo.
The `dbfv` and `dckks` packages import the `drlwe` and provide decorators over its types.

## MHE-MPC Protocol Overview

The MHE-MPC protocol implemented in Lattigo is based on the constructions described in ["Multiparty Homomorphic Encryption from Ring-Learning-with-Errors"](https://eprint.iacr.org/2020/304.pdf) by Mouchet et al. (2021), which is a RLWE instantiation of the MPC protocol described in ["Multiparty computation with low communication, computation and interaction via threshold FHE"](https://eprint.iacr.org/2011/613.pdf) by Asharov et al. (2012).

The protocol permits to a group of parties to compute a joint function of their private inputs under encryption and to release the result in plaintext to an arbitrary receiver.
The protocol provides security against *passive* attackers that can be the parties themselves and can operate in various system models:

**Anytrust vs Threshold Access-structure**. In the traditional MPC setting, the parties encrypt their data such that no coalition of at most N-1 parties can decrypt which preserves confidentiality of the inputs as long as at least one party is honest. However, this might be too restrictive in some setting where at least _t_ parties can be assumed honest. The protocol support _t-out-of-N_ threshold types of access-structure for the latter case.

**Peer-to-peer vs Cloud-assisted models**. The parties can execute the protocol in a peer-to-peer way or be assisted by a third-party server (which is also considered a passive adversary).

**Internal vs External receivers**. _Receivers_ are the intended recipients of the computation and can be either _internal_ or _external_ depending or whether they are or are not input parties themselves, respectively. The MHE-MPC protocol provides two modes for its Output procedure, one for each type of receiver.

An execution of the complete MHE-based MPC protocol has the following structure (the details of each protocols are provided later):
1. Setup  
  1. Secret Keys Generation        
  1. [OPTIONAL] Threshold Secret-Key Generation
  1. Collective Public Encryption-Key Generation
  1. [OPTIONAL] Collective Public Relinearization-Key Generation  
  1. [OPTIONAL] Collective Public Rotation-Key Generation
2. Input Phase (Encryption)
3. Evaluation Phase
4. Output phase  (Decryption)
  1. [OPTIONAL] Threshold Combine Phase  
  1. Collective Key-Switching
  1. Local Decryption
___
## MHE-MPC Protocol Steps Description
This section provides a description for each sub-protocol composing the MHE-MPC protocol. 
The system model is abstracted by considering that the parties have access to a common public authenticated channel.
In the cloud-assisted setting, this public channel can be the helper cloud-server.
In the peer-to-peer setting , it could be a broadcast channel or a topological one.

### Setup
In this phase, the parties generate the various keys that are required by the Input, Evaluation and Output phases.

#### Secret Keys Generation
The parties generate their individual secret-keys locally using a `rlwe.KeyGenerator`; this provides them with a `rlwe.SecretKey` type.
See [rlwe/keygen.go](../rlwe/keygen.go) further information on key-generation.
The _ideal secret-key_ is defined as the sum of all secret-keys.

This secret-key enforces an _N-out-N_ access structure which requires all the parties to collaborate in a ciphertext decryption (hence, tolerates N-1 dishonest parties).

#### [OPTIONAL] Threshold Secret-Key Generation
For settings where an _N-out-N_ access structure is too restrictive (e.g., from an availability point of view), an optional Threshold Secret-Key Generation Protocol can be performed to enable _t-out-of-N_ access-structures (hence tolerating t-1 dishonest parties).
The idea of this protocol is to apply Shamir Secret Sharing to the _ideal secret-key_ in such a way that any group of _t_ parties can reconstruct it.
This is achieved by a single-round protocol where each party applies Shamir Secret-Sharing to its own share of the _ideal secret-key_. The sharing is done in the natural domain of these share.

We assume that each party is associated with a `drlwe.ShamirPublicKey` that is known to the other parties.

This protocol is implemented through the `drlwe.ThresholdizerProtocol` interface and its steps are as follows:
- Each party generates a `drlwe.ShamirPolynomial` by using the `drlwe.ThresholdizerProtocol.GenShamirPolynomial` method, then generates a share of type `drlwe.ShamirSecretShare` for each of the other parties' `ShamirPublicKey` by using the `drlwe.ThresholdizerProtocol.GenShamirSecretShare`.
- Every party use the IDs to compute each party's "threshold public key", which is a point in the ring given at the ThresholdizerProtocol initiation, that has a type `ThreshPublicKey`.
- Each party sends the respective `ShamirSecretShare` to each of the other parties over a private channel. 
- Each party aggregates the all `ShamirSecretShare`s it receives using the `drlwe.ThresholdizerProtocol.AggregateShares` method.

Each party stores its aggregated `ShamirSecretShare` for later use in the Combining step of the Output Phase.

#### Public Key Generation
The parties execute the collective public encryption-key generation protocol to obtain an encryption-key for the _ideal secret-key_.

The protocol is implemented through the `drlwe.CollectivePublicKeyGenerator` interface and its steps are as follows:
- Each party generates a share (`drlwe.CKGShare`) from their secret-key, by using the `drlwe.CollectivePublicKeyGenerator.GenShare` method.
- Each party discloses its share over the public channel. The shares are aggregated with the `drlwe.CollectivePublicKeyGenerator.AggregateShare` method.
- Each party can derive the public encryption-key (`rlwe.PublicKey`) by using the `drlwe.CollectivePublicKeyGenerator.GenPublicKey` method.

Note that the `GenShare` requires a polynomial that is uniformly random and common to all parties (`crs`).

#### [OPTIONAL] Relinearization Key Generation
This protocol provides the parties with a public relinearization-key (`rlwe.RelinearizationKey`) for the _ideal secret-key_. This public-key enables compact multiplications in RLWE schemes. If the circuit contains no multiplication (e.g. a purely additive circuit), no relinearization key is needed. It has two rounds.

The protocol is implemented through the  `drlwe.RelinearizationKeyGenerator` interface and its steps are as follows:
- Each party generates a share (`drlwe.RGKShare`) for the first protocol round by using the `drlwe.RelinearizationKeyGenerator.GenShareRoundOne` method. This method also provides the party with an ephemeral secret-key (`rlwe.SecretKey`) that is required for the second round.
- Each party discloses its share for the first round over the public channel. The shares are aggregated with the `drlwe.RelinearizationKeyGenerator.AggregateShare` method.
- Each party generates a share (also a `drlwe.RGKShare`) for the second protocol round by using the `drlwe.RelinearizationKeyGenerator.GenShareRoundTwo` method.
- Each party discloses its share for the second round over the public channel. The shares are aggregated with the `drlwe.RelinearizationKeyGenerator.AggregateShare` method.
- Each party can derive the public relinearization-key (`rlwe.RelinearizationKey`) by using the `drlwe.RelinearizationKeyGenerator.GenRelinearizationKey` method.

Note that, similar to the CKG protocol, the RKG protocol requires a matrix of polynomials that are uniformly random and common to all parties (`crs`).

#### [OPTIONAL] Rotation Key Generation
This protocol provides the parties with a public rotation-key (stored as `rlwe.SwitchingKey` types) for the _ideal secret-key_. One rotation-key enables one specific rotation on the ciphertexts' slots. The protocol can be repeated to generate the keys for multiple rotations.

The protocol is implemented through the  `drlwe.RotationKeyGenerator` interface and its steps are as follows:
- Each party generates a share (`drlwe.RTGShare`) by using `drlwe.RotationKeyGenerator.GenShare`. 
- Each party discloses its `drlwe.RTGShare` over the public channel. The shares are aggregated with the `drlwe.RotationKeyGenerator.AggregateShare` method.
- Each party can derive the public rotation-key (`rlwe.SwitchingKey`) from the final `RTGShare` by using the `drlwe.RotationKeyGenerator.AggregateShare` method.

Note that, similar to the CKG and RKG protocols, the RKG protocol requires a matrix of polynomials that are uniformly random and common to all parties (`crs`).

---
### Input Phase (Encryption Phase)
The parties provide their inputs for the computation during the Input Phase.
They use the collective encryption-key generated during the Setup Phase to encrypt their inputs, and disclose them on the public channel.
Since the collective encryption-key is a valid RLWE public encryption-key, it can be used directly with the single-party scheme.
Hence, the parties use the `Encoder` and `Encryptor` interfaces of the desired encryption scheme (see [bfv.Encoder](../bfv/encoder.go), [bfv.Encryptor](../bfv/encryptor.go), [ckks.Encoder](../ckks/encoder.go) and [ckks.Encryptor](../ckks/encryptor.go)).

---
### Evaluation phase
The computation of the desired function is performed homomorphically during the Evaluation Phase.
The step can be performed by the parties themselves or can be outsourced to a cloud-server. 
Since the ciphertexts in the multiparty schemes are valid ciphertexts for the single-party ones, the homomorphic operation of the latter can be used directly (see [bfv.Evaluator](../bfv/evaluator.go) and [ckks.Evaluator](../ckks/evaluator.go)).

---
### Output phase
The receiver(s) obtain their outputs through the final Output Phase, which aim is to decrypt the ciphertexts resulting from the Evaluation Phase.
It is a two-steps process with an optional pre-processing step in case of the t-out-of-N access-structure.
In the first step, Collective Key-Switching the parties re-encrypt the desired ciphertext under the receiver's secret-key.
The second step is the local decryption of this re-encrypted ciphertext by the receiver.

#### Combining Step (Optional)
If not all parties are online during the Output phase, but a Threshold Key-Generation step was performed, the parties can perform a Combining step to compute a new sharing of the _ideal secret-key_ with only _t_ parts. This step is local to all parties and does not require interaction once the set of online parties is known.

This step is implemented by the `drlwe.Combiner` interface, which simply provides `GenAdditiveShare` method that computes the party's new share of the _ideal secret-key_ given the list of `ShamirPublicKey` of the active parties and the party's `ShamirSecretShare`. The result is stored in a `rlwe.SecretKey` that can be used in place of the _N-out-of-N_ one during the next phase performed among _t_ parties.

#### Collective Key-Switching
The parties perform a re-encryption of the desired ciphertext(s) from being encrypted under the _ideal secret-key_ to being encrypted under the receiver's secret-key.
There are two instantiation of the Collective Key-Switching protocol:
- Collective Key-Switching (CKS), implemented as the `drlwe.KeySwitchingProtocol` interface, enables the parties to switch from their _ideal secret-key_ _s_ to another _ideal secret-key_ _s'_ when s' is collectively known to the parties. In the case where _s' = 0_, this is equivalent to a decryption protocol which can be used when the receiver is one of the input-parties. 
- Collective Public-Key Switching (PCKS), implemented as the `drlwe.PublicKeySwitchingProtocol` interface, enables parties to switch from their _ideal secret-key_ _s_ to an arbitrary key _s'_ when provided with a public encryption-key for _s'_. Hence, this enables key-switching to a secret-key that is not known to the input parties, which further enables external receivers.

While both protocol variants have slightly different local operation, their steps are the same and as follows:
- Each party generates a share (of type `drlwe.CKSShare` or `drlwe.PCKSShare`) with the `drlwe.(Public)KeySwitchingProtocol.GenShare` method. This requires its own secret-key (a `rlwe.SecretKey`) as well as the destination key: its own share of the destination key (a `rlwe.SecretKey`) in CKS or the destination public-key (a `rlwe.PublicKey`) in PCKS.
- Each party discloses its `drlwe.CKSShare` over the public channel. The shares are aggregated with the `drlwe.(Public)KeySwitchingProtocol.AggregateShares` method.
- From the aggregated `drlwe.CKSShare`, any party can derive the ciphertext re-encrypted under _s'_ by using the `drlwe.(Public)KeySwitchingProtocol.KeySwitch` method.

#### Decryption
Once the receivers have obtained the ciphertext re-encrypted under their respective keys, they can use the usual decryption algorithm of the single-party scheme to obtain the plaintext result (see [bfv.Decryptor](../bfv/decryptor.go) and [ckks.Decryptor](../ckks/decryptor.go)).

