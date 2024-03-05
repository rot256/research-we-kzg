# KZG Witness Encryption

This is the implementation artifacts for the paper:

[Extractable Witness Encryption for KZG Commitments and Efficient Laconic OT](https://eprint.iacr.org/2024/264)

**Abstract:**
We present a concretely efficient and simple extractable witness encryption scheme for KZG polynomial commitments.
It allows to encrypt a message towards a triple $(\mathsf{com}, \alpha, \beta)$, where $\mathsf{com}$ is a KZG commitment for some polynomial $f(X)$.
Anyone with an opening for the commitment attesting $f(\alpha) = \beta$ can decrypt, but without knowledge of a valid opening the message is computationally hidden.

Our construction is simple and highly efficient. The ciphertext is only a single group element. Encryption and decryption both require a single pairing evaluation and a constant number of group operations.

Using our witness encryption scheme, we construct a simple and highly efficient laconic OT protocol, which significantly outperforms the state of the art in most important metrics.

## Reproduction

Benchmarks can be reproduced by simply running `cargo bench`.
