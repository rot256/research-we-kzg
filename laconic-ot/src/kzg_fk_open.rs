use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup};
use ark_poly::EvaluationDomain;
use ark_std::Zero;
use std::ops::Mul;

use crate::kzg_types::{CommitmentKey, State};

// this module allows to compute all openings in a
// fast amortized way following the FK technique:
// https://eprint.iacr.org/2023/033.pdf

/// Compute the vector y = DFT(hat_s) used for fast amortized
/// computation of all KZG openings. It is independent of the
/// committed vector and can therefore be computed once from
/// the public parameters / commitment key
pub fn precompute_y<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    powers: &[E::G1Affine],
    domain: &D,
) -> Vec<E::G1Affine> {
    // Preparation: We need to do FFTs of twice the size.
    // for that, we make use of a evaluation domain of twice the size
    // we also let d be the degree of f
    let d = domain.size() - 1;
    let domain2 = D::new(2 * domain.size()).unwrap();
    // hat_s = [powers[d-1],...,powers[0], d+2 neutral elements]
    let mut hat_s = Vec::with_capacity(2 * d + 2);
    for i in (0..=d - 1).rev() {
        hat_s.push(powers[i].into_group());
    }
    for _ in 0..d + 2 {
        hat_s.push(E::G1::zero());
    }
    let y = domain2.fft(&hat_s);
    E::G1::normalize_batch(&y)
}

/// FK technique to compute openings in a *non-hiding* way
/// evals contains the domain.size() many evaluations
/// of the polynomial over the evaluation domain
pub fn all_openings_single<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    y: &[E::G1Affine],
    domain: &D,
    evals: &[E::ScalarField],
) -> Vec<E::G1> {
    // compute the base polynomial h
    let coeffs = domain.ifft(&evals);
    let mut h = base_poly::<E, D>(y, domain, &coeffs);

    // evaluate h in the exponent using FFT
    // the evaluations are the openings
    domain.fft_in_place(&mut h);
    h
}

/// compute the polynomial h (in exponent) from the paper (see Proposition 1)
/// The polynomial f is given by domain.size() many coefficients, and we have
/// powers[i] = g1^{alpha^i}
/// The ith KZG opening is h(domain.element(i)). Hence, one we have h, we can
/// compute all openings efficiently using a single FFT in the exponent
fn base_poly<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    y: &[E::G1Affine],
    domain: &D,
    coeffs: &[E::ScalarField],
) -> Vec<E::G1> {
    // we follow the modifications as in the implementation of caulk
    // https://github.com/caulk-crypto/caulk/blob/main/src/dft.rs#L17

    // Preparation: We need to do FFTs of twice the size.
    // for that, we make use of a evaluation domain of twice the size
    // we also let d be the degree of f
    let d = domain.size() - 1;
    let domain2 = D::new(2 * domain.size()).unwrap();

    // Step 1: y = DFT(hat_s) has already been precomputed
    // Step 2: v = DFT(hat_c), where

    // hat_c = [coeffs[d], d zeros, coeffs[d], coeffs[0],...,coeffs[d-1]]
    let mut hat_c = Vec::with_capacity(2 * d + 2);
    hat_c.push(coeffs[d]);
    for _ in 0..d {
        hat_c.push(E::ScalarField::zero());
    }
    hat_c.push(coeffs[d]);
    for i in 0..d {
        hat_c.push(coeffs[i]);
    }

    // let v = domain2.fft(&hat_c);
    domain2.fft_in_place(&mut hat_c);
    let v = hat_c;

    // Step 3: u = comp.-wise prod. of y and v
    let mut u: Vec<E::G1> = Vec::with_capacity(2 * d);
    for i in 0..2 * d + 2 {
        u.push(y[i].mul(v[i]));
    }

    // Step 4: hat_h = iDFT(u)
    //let hat_h = domain2.ifft(&u);
    domain2.ifft_in_place(&mut u);
    let hat_h = u;
    let h = hat_h[0..d].to_vec();
    h
}

#[cfg(test)]
mod tests {

    use std::vec;

    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_ec::{CurveGroup, VariableBaseMSM};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::EvaluationDomain;
    use ark_poly::{DenseUVPolynomial, Radix2EvaluationDomain};
    use ark_std::One;
    use ark_std::UniformRand;

    use crate::kzg_types::VcKZG;

    use super::{all_openings_single, base_poly};

    type F = <Bls12_381 as Pairing>::ScalarField;
    type D = Radix2EvaluationDomain<F>;

    /// test function base_polynomial
    #[test]
    fn test_base_poly() {
        let mut rng = ark_std::rand::thread_rng();
        let degree = 15;
        let runs = 10;

        // generate some parameters
        let ck = VcKZG::<Bls12_381, D>::setup(&mut rng, degree - 1).unwrap();

        for _ in 0..runs {
            // sample random polynomial f and its evaluations
            let f = DensePolynomial::rand(ck.domain.size() - 1, &mut rng);
            // compute the expected coefficients of h
            // in the exponent (expensive version)
            let mut naive = Vec::new();
            for i in 1..=degree {
                // according to paper:
                // h_i = f[d]u[d-i] + f[d-1]u[d-i-1] + ... + f[i+1]u[1] + f[i]u[0],
                // where u[j] = g1^{secret^j} and d = degree
                // note that this is an MSM of f[i..=d] and u[0..=d-i]
                let hi = <<Bls12_381 as Pairing>::G1 as VariableBaseMSM>::msm(
                    &ck.u[0..=(degree - i)],
                    &f.coeffs[i..=degree],
                )
                .unwrap()
                .into_affine();
                naive.push(hi);
            }
            // compute h using the function we want to test
            let h = base_poly::<Bls12_381, D>(&ck.y, &ck.domain, &f.coeffs);
            // check that they are indeed equal
            for i in 0..=degree - 1 {
                assert_eq!(naive[i], h[i].into_affine());
            }
        }
    }

    /// test function all_openings_single
    #[test]
    fn test_all_openings_single() {
        let mut rng = ark_std::rand::thread_rng();
        let degree = 15;
        let runs = 10;

        // generate some parameters
        let ck = VcKZG::<Bls12_381, D>::setup(&mut rng, degree - 1).unwrap();

        for _ in 0..runs {
            // generate random polynomial and its evaluations
            let f = DensePolynomial::rand(ck.domain.size() - 1, &mut rng);
            let evals = ck.domain.fft(&f.coeffs);
            // precompute the openings naively using long division (very slow)
            let mut naive: Vec<<Bls12_381 as Pairing>::G1Affine> = Vec::new();
            for i in 0..ck.domain.size() {
                // witness poly using long division
                let z = ck.domain.element(i);
                let fshift = &f - &DensePolynomial::from_coefficients_vec(vec![evals[i]]);
                let div = DensePolynomial::from_coefficients_vec(vec![-z, F::one()]);
                let witness_poly = &fshift / &div;
                // commit to witness poly at alpha
                let c = <<Bls12_381 as Pairing>::G1 as VariableBaseMSM>::msm(
                    &ck.u[0..ck.domain.size() - 1],
                    &witness_poly.coeffs,
                )
                .unwrap();
                naive.push(c.into_affine());
            }
            // precompute the openings using the function we want to test
            let fk: Vec<<Bls12_381 as Pairing>::G1> =
                all_openings_single::<Bls12_381, D>(&ck.y, &ck.domain, &evals);
            // compare the results
            for i in 0..ck.domain.size() {
                assert_eq!(naive[i], fk[i].into_affine());
            }
        }
    }

    // test the public function all_openings
    #[test]
    fn test_all_openings() {
        let mut rng = ark_std::rand::thread_rng();
        let degree = 15;
        let runs = 3;

        for _ in 0..runs {
            // generate some parameters
            let ck = VcKZG::<Bls12_381, D>::setup(&mut rng, degree - 1).unwrap();

            // commit to something
            let m = (0..degree - 1).map(|_| F::rand(&mut rng)).collect();
            let (_com, mut st) = VcKZG::<Bls12_381, D>::commit(&mut rng, &ck, &m);

            // compute all the openings freshly
            let mut openings = Vec::new();
            for i in 0..ck.message_length {
                let op = VcKZG::<Bls12_381, D>::open(&ck, &st, i as u32).unwrap();
                openings.push(op.v);
            }

            // check that all openings are the same
            let precomputed = st.precomputed_v.unwrap();
            for i in 0..ck.message_length {
                assert_eq!(precomputed[i], openings[i]);
            }
        }
    }
}
