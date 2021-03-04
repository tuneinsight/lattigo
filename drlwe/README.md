# DRLWE
The DRLWE package is a collection of generalized protocols based on the ring-learning-with-error problem, in order to evaluate a function on distributed inputs.

## Aim
The aim of the package is to provide generic interfaces for different steps of a distributed homomorphic encryption protocol, that can then be implemented in one of the schemes supported by Lattigo (DBFV of DCKKS).

## General protocols
An execution of the complete scheme has the following protocol structure, with some protocols being optional. The details of each protocol are provided later.
1. Setup  
  1. Secret Keys Generation        
  1. Public Key Generation
  1. [OPTIONAL] Threshold share phase  
  1. [OPTIONAL] Relinearization Key Generation  
  1. [OPTIONAL] Rotation Key Generation
2. Input Phase (Encryption Phase)
3. Evaluation Phase
4. Out phase  
  1. [OPTIONAL] Threshold Combine Phase  
  1. (Public) Key Switching  
  1. Decryption
___
## Protocols explained
### Setup
#### Secret Keys Generation
At the very start of any execution of the scheme, each party must receive or compute a secret key of type `rlwe.SecretKey`. This is to be handled by the user using for example a `KeyGenerator`.

See [bfv.keygen](../bfv/keygen.go) and [ckks.keygen](../ckks/keygen.go) for further information.

#### Public Key Generation
The goal of this protocol is to generate a common public key from the secret of all parties.

This protocol requires a Common Reference String (CRS), which is a string known by every party. The protocol steps are as follows:
- Every party generates a share from their secret s (of type `CKGShare`) of the form ((crs*s + e)), where e is sampled from a gaussian distribution with low magnitude.
- Every party broadcasts its share. Because of the assumptions on the ring-learning-with-error, an adversary cannot infer s from the value of the share.
- Every party aggregates (by adding) all the share they have received. The final aggregate is the same for every party, and is the public key, of type `rlwe.PublicKey`.

The methods offered by a Common Key Generation Protocol are:
#### CKGProtocol
- `AllocateShares`, to create an empty `CKGShare`.                
- `GenShare`, to generate a share as described before.
- `AggregateShare`, to aggregate 2 `CKGShare` into 1 by adding them.
- `GenPublicKey`, to convert the final `CKGShare` aggregate into a `rlwe.PublicKey`

#### [OPTIONAL] Threshold Sharing Phase
If not all parties are guaranteed to be operational at decryption time, the user has the possibility to use a t-out-of-N setup with the help of a thresholdizing protocol, and later a combining protocol.

This protocol requires every party to have a unique, public identifier (ID). the generation and distribution of these IDs is left to the user.
The sharing phase of the threshold protocol has the following steps :
- Every party samples a random polynomial whose constant term is the value of their secret key. Such a polynomial is called a "share generating polynomial" and has a type `ShareGenPoly`.
- Every party use the IDs to compute each party's "threshold public key", which is a point in the ring given at the ThresholdizerProtocol initiation, that has a type `ThreshPublicKey`.
- Every party evaluates their share generating polynomial at other parties' threshold public keys. The result of this evaluation is called a "threshold secret share" (and has a type `ThreshSecretShare`), and must be sent to the corresponding key holders.
- Every party aggregates all the shares it receives into what we call a "threshold secret key", which has type `rlwe.SecretKey`.

Each party should store their threshold secret key as well as all parties' threshold public keys until the combining phase.

The methods offered by a Threshold Sharing Protocol (also called Thresholdizer Protocol) are detailed in the following:
###### ThresholdizerProtocol
An implementation of a `ThresholdizerProtocol` offers the following methods :
- `Allocate*`, methods to create empty instances of objects.
- `GenKeyFromID` that takes a `PublicID` and returns a `ThreshPublicKey`, such that the result is deterministic and consistent among different Go routines.
- `InitShareGenPoly`, that creates a slice of coefficients such that the first is the party's secret share, and the other are uniformly random.
- `GenShareForParty`, that use a `ShareGenPoly` to generate a share for a given party by evaluating the `ShareGenPoly`'s polynomial at the point corresponding to the given party's `ThreshPublicKey`
- `GenThreshSecretKey`, that generates a `rlwe.SecretKey` from an aggregation of `ThresholdSecretShare`
- `AggregateShares`, to aggregate two `ThreshSecretShare` into one.

#### [OPTIONAL] Relinearization Key Generation
The aim of this protocol is to provide the scheme with a relinearization key (of type `rlwe.RelinearizationKey`), in order to integrate multiplications. If the circuit contains no multiplication (e.g. a purely additive circuit), no relinearization key is needed.

Similar to the CKG protocol, the RKG protocol requires a Common Reference Polynomial (CRP). Also, each party needs an additional, ephemeral secret key (of type `rlwe.SecretKey`)
The Relinearization Key Generation Protocol has the following steps:
- First round of share generation : each party generates an ephemereal key and uses it to pseudo-encrypt its secret key share. This generates a share `RGKShare`, that is broadcast.
- Each party aggregates the shares it receives from the first round into a combined share.
- Second round of share generation: each party uses their ephemeral and secret keys as well as the combined share from round one to generate a second share, of type `RKGShare`.
- Each party aggregates the shares it receives from the first round into a second combined share.
- Each party uses the first and second combined share to generate the `RelinearizationKey`.

The methods offered by a RKG Protocol are the following:

##### RelinearizationKeyGenerator
- `AllocateShares()`, a method to create empty instances of an ephemeral `SecretKey` and two `RKGShare`.
- `GenShareRoundOne`, that takes the party's `SecretKey` and the CRP to produce the first round `RKGShare` as well as an ephemeral `SecretKey`.
- `GenShareRoundTwo` uses the party's `SecretKey` as well as its ehpemeral `SecretKey`, the CRP and the aggregation of all received `RKGShare` of round 1 to generate another `RKGShare`.
- `AggregateShares` aggregates two `RKGShare` into one.
- `GenRelinearizationKey` takes the aggregation of both rounds to generate a `RelinearizationKey`

#### [OPTIONAL] Rotation Key Generation
The aim of this protocol is to provide the scheme with a rotation key (of type `rlwe.SwitchingKey`), in order to enable the parties to perform rotations on ciphertexts, to evaluate non-linear functions.

As for the RKG Protocol, the Rotation Key Generation (RTG) Protocol needs a Common Reference Polynomial (CRP), that may differ from the one used in an eventual prior RKG.
A RTG Protocol is very similar to a CKG Protocol and has the following steps:

- Each party generates a share (of type `RTGShare`) from its `SecretKey`, the CRP and a Galois element (?).
- Each party broadcasts its `RTGShare`.
- Each party aggregates every `RTGShare` it received into a final combined share, that is the same for every party. This final share is converted to a `RotationKey`.

A Rotation Key Generation Protocol offers the following methods :

##### RotationKeyGenerator
- `AllocateShares`, that creates an empty instance of `RTGShare`
- `GenShare`, that generates a `RTGShare` from a party's `SecretKey`, a Galois element (of type `uint64`), and the CRP.
- `Aggregate`, that aggregates two `RTGShare` into one.
- `GenRotationKey`, that generates a rotation key of type `rlwe.SwitchingKey` from the aggregated `RTGShare` and the CRP.

---
### Input Phase (Encryption Phase)
The aim of this protocol is for every party to encrypt their input using the public key generated at step 3.

To proceed to this step, the user might want to use an `Encoder` and/or an `Encryptor`, that are not part of this package.

See [bfv.encoder](../bfv/encoder.go), [bfv.encryptor](../bfv/encryptor.go), [ckks.encoder](../ckks/encoder.go) and [ckks.encryptor](../ckks/encryptor.go) for further information.

---
### Evaluation phase
The aim of this protocol is to produce a ciphertext that is the result of the
desired function on the parties' ciphertexts. To proceed to this step, the user might want to use an `Evaluator`, that is not part of this package.

See [bfv.evaluator](../bfv/evaluator.go) and [ckks.evaluator](../ckks/evaluator.go) for further information.

---
### Out phase
#### Threshold Combining Phase (Optional)
This protocol is mandatory if the Threshold Share phase has been executed.

It has the following steps:
- The cohort has to determine which parties are active.
- **Active parties** compute the Lagrange interpolation of their threshold public key with all other **active players'** public keys, with a coefficient corresponding to their threshold secret key. The result is a share that must be used **instead** of the party's secret key share for the key switching.

Note that this step does not require any communication.

The methods offered by a Combiner Protocol are :
##### CombinerProtocol
- `GenFinalShare` : from the list of active players' `ThresholdPublicKey`, the party's `ThresholdPublicKey` as well as the party's threshold `SecretKey` that was generated during threshold sharing, generate a share that can be used by this party for public switching instead of the standard secret key share of this party.
- `Equal` : compares two `ThresholdPublicKey` for equality. This method ought to be efficient as it will be called many times during Lagrange Interpolation.

#### (Public - ) Collective Key switching
The aim of this protocol is, from a ciphertext encrypted with the cohort's public key, to produce a ciphertext encrypted with another (arbitrary) public key, such that the decryption of the two ciphertexts using their respective secret keys is equal.

The difference between a Collective Key Switching (CKS) Protocol and a Public Collective Key Switching (PCKS) Protocol is that the first assumes the cohort has collective knowledge of the output secret key, whereas the second only relies on a given output public key.

The steps of both versions are pretty much the same and are as follows:
- Each party generates a share (of type `CKSShare` (CKS) or `PCKSShare` (PCKS) from its `SecretKey`, the input ciphertext and the output secret key (KCS) or the output public key (PCKS).
- Each party aggregates every `(P)CKSShare` it receives into a final combined share, that is the same for every party.
- From the combined share, and the input ciphertext, each party perform the switch and outputs a ciphertext that can be decrypted using the output secret key (CKS) or the secret key corresponding to the output public key (PCKS).

A Key swithing protocol offers the following methods:
- `AllocateShare`, to create an empty instance of `(P)CKSShare`
- `GenShare`, to generate a `(P)CKSShare` from a `SecretKey`, the input `Ciphertext` and either the output `SecretKey` (CKS) or the output `PublicKey` (PCKS)
- `AggregateShares`, to combine 2 `(P)CKSShare` into 1.
- `KeySwitch`, to perform the key switching from the aggregated `(P)CKSShare` and the input `Ciphertext`, producing the output `Ciphertext`.

#### Decryption
The aim of this protocol is to get a plaintext from a ciphertext and a secret key.

The user might want to use a `Decryptor`, which is not provided in this package.

See [bfv.decryptor](../bfv/decryptor.go) and [ckks.decryptor](../ckks/decryptor.go) for further information.
