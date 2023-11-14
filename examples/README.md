# Single Party Examples

## Applications

Application examples are examples showcasing specific capabilities of the library or scaled down real world scenarios.

### Binary

- `bin_blind_rotations`: an example showcasing the evaluation of the sign function using blind rotations on RLWE ciphertexts.

### Integers

- `int_ride_hailing`: an example on privacy preserving ride hailing.
- `int_vectorized_OLE`: an example on vectorized oblivious linear evaluation using RLWE trapdoor.

### Reals/Complexes

- `reals_bootstrapping`: a series of example showcasing the capabilities of the bootstrapping for fixed point arithmetic.
  - `basics`: an example showcasing the basic capabilities of the bootstrapping.
  - `high_precision`: an example showcasing high precision bootstrapping.
  - `slim`: an example showcasing slim bootstrapping, i.e. re-ordering the steps of the bootstrapping.

- `reals_scheme_switching`: an example showcasing scheme switching between `hefloat` and `hebin` to complement fixed-point arithmetic with lookup tables.
- `reals_sigmoid_chebyshev`: an example showcasing polynomial evaluation of a Chebyshev approximation of the sigmoid.
- `reals_sigmoid_minimax`: an example showcasing polynomial evaluation of a minimax approximation of the sigmoid.
- `reals_vectorized_polynomial_evaluation`: an example showcasing vectorized polynomial evaluation, i.e. evaluating different polynomials in parallel on specific slots.

## Templates

Templates are files containing the basic instantiation, i.e. parameters, key-generation, encoding, encryption and decryption.

- `reals`: a template for `hefloat`.
- `int`: a template for `heint`.

## Tutorials

Tutorials are examples showcasing the basic capabilities of the library.

- `reals`: a tutorial on all the basic capabilities of the package `hefloat`.

# Multi Party Examples

 - `int_pir`: an example showcasing multi-party private information retrieval.
 - `int_psi`: an example showcasing multi-party private set intersection.
 - `thresh_eval_key_gen`: an example showcasing multi-party threshold key-generation.