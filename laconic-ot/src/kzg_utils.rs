use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::{batch_inversion, Field};
use ark_poly::{EvaluationDomain, Polynomial};
use ark_serialize::CanonicalSerialize;
use ark_std::{One, Zero};
use std::ops::Mul;

use crate::kzg_types::{Commitment, CommitmentKey, Opening};

// This module contains helper functions for the Simulation Extractable KZG Vector commitment

/// Standard KZG verification. Verifies that f(z) = y
#[inline]
pub fn plain_kzg_verify<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    ck: &CommitmentKey<E, D>,
    com_kzg: &E::G1Affine,
    z: E::ScalarField,
    y: E::ScalarField,
    tau: &Opening<E>,
) -> bool {
    // check e(com*g1^{-y}*h^{-hat_y},g2) == e(v,r*g2^{-z})
    let mut lhs_left = com_kzg.into_group();
    lhs_left -= ck.u[0].mul(y);
    let rhs_right = ck.r.into_group() - ck.g2.mul(z);
    // Naive Implementation:
    //  let lhs = E::pairing(lhs_left, ck.g2);
    //  let rhs = E::pairing(tau.v, rhs_right);
    //  lhs == rhs
    // We can do it slightly faster:
    let left = vec![E::G1Prepared::from(-lhs_left), E::G1Prepared::from(tau.v)];
    let right = vec![E::G2Prepared::from(ck.g2), E::G2Prepared::from(rhs_right)];
    let q = E::multi_pairing(left, right);
    q.is_zero()
}

/// Standard KZG verification. Verifies that f(z) = y,
/// but assumes that z = domain.element(i)
#[inline]
pub fn plain_kzg_verify_inside<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    ck: &CommitmentKey<E, D>,
    i: usize,
    com_kzg: &E::G1Affine,
    y: E::ScalarField,
    tau: &Opening<E>,
) -> bool {
    // check e(com*g1^{-y}*h^{-hat_y0},g2) == e(v0,r*g2^{-z0})
    let mut lhs_left = com_kzg.into_group();
    lhs_left -= ck.u[0].mul(y);
    // Naive Implementation:
    //  let lhs = E::pairing(lhs_left, ck.g2);
    //  let rhs = E::pairing(tau.v, ck.d[i]);
    //  lhs == rhs
    // We can do it slightly faster:
    let left = vec![E::G1Prepared::from(-lhs_left), E::G1Prepared::from(tau.v)];
    let right = vec![E::G2Prepared::from(ck.g2), E::G2Prepared::from(ck.d[i])];
    let q = E::multi_pairing(left, right);
    q.is_zero()
}

/// Compute a KZG commitment for the given vector of evaluations
#[inline]
pub fn plain_kzg_com<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    ck: &CommitmentKey<E, D>,
    evals: &[E::ScalarField],
) -> E::G1Affine {
    assert_eq!(evals.len(), ck.lagranges.len());
    let c = <E::G1 as VariableBaseMSM>::msm(&ck.lagranges, evals).unwrap();
    c.into_affine()
}

/// Check if the given element is in the evaluation domain
/// and if so, return the index of it. Otherwise, return None
#[inline]
pub fn find_in_domain<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    domain: &D,
    z: E::ScalarField,
) -> Option<usize> {
    if domain.vanishing_polynomial().evaluate(&z) == E::ScalarField::zero() {
        // this should happen with negl prob.
        let mut y = None;
        for i in 0..domain.size() {
            if z == domain.element(i) {
                y = Some(i);
            }
        }
        y
    } else {
        None
    }
}

/// Compute the evaluation form of the KZG witness polynomial
/// psi = (f - f(w_i)) / (X - w_i) when f is given in evaluation form
/// Note: This assumes that w_i is the ith element of the domain
/// The evaluation form is appended to the given vector witn_evals
/// that is, the jth pushed element is psi(w_j)
#[inline]
pub fn witness_evals_inside<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    domain: &D,
    evals: &[E::ScalarField],
    i: usize,
    witn_evals: &mut Vec<E::ScalarField>,
) {
    // need that later for index calculation
    let oldsize = witn_evals.len();

    // let x_j denote the elements of the evaluation domain
    // then for each j != i, we can compute the evaluation
    // as in witness_evals_outside,
    // namely, as (evals[j] - evals[i]) / (x_j-x_i)
    // for that, we first compute all
    // denominators we need using batch inversion
    let fxi = evals[i];
    let xi = domain.element(i);
    let mut nums = Vec::new();
    let mut denoms = Vec::new();
    for j in 0..domain.size() {
        // f(x_j) - f(x_i)
        nums.push(evals[j] - fxi);
        // x_j-x_i
        denoms.push(domain.element(j) - xi);
    }
    // now, denoms[i] = 0. So let's set it to 1
    // to make batch inversion possible
    denoms[i] = E::ScalarField::one();
    batch_inversion(&mut denoms);
    for j in 0..domain.size() {
        witn_evals.push(nums[j] * denoms[j]);
    }
    // now witn_evals is correctly computed for all j!=i.
    // whats left is to compute the ith evaluation properly
    // https://dankradfeist.de/ethereum/2021/06/18/pcs-multiproofs.html
    witn_evals[oldsize + i] = {
        let mut sum = E::ScalarField::zero();
        for j in 0..domain.size() {
            if j == i {
                continue;
            }
            let mut term = nums[j] * (-denoms[j]);
            let d = domain.size();
            let exponent = (j as isize) - (i as isize);
            let exponent = ((exponent + d as isize) as usize) % d;
            term *= domain.element(exponent);
            sum += term;
        }
        sum
    };
}

/// computes the vector of all 1/(domain[i]-z)
/// Assumes that z is not in domain
#[inline]
pub fn inv_diffs<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    domain: &D,
    z: E::ScalarField,
) -> Vec<E::ScalarField> {
    // we use batch inversion for the denominators
    let mut inv_diffs = Vec::with_capacity(domain.size());
    for i in 0..domain.size() {
        inv_diffs.push(domain.element(i) - z);
    }
    batch_inversion(&mut inv_diffs);
    inv_diffs
}

/// Compute the evaluation form of the KZG witness polynomial
/// Evaluation form will be pushed to the vector witn_evals
/// ith pushed element is (f(domain[i]) - f(z)) / (domain[i] - z)
/// where i ranges from 0 to domain.size()
/// Assumes inv_diffs[i] = 1/(domain[i]-z) for i in 0..domain.size()
#[inline]
pub fn witness_evals_outside<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    domain: &D,
    evals: &[E::ScalarField],
    fz: E::ScalarField,
    inv_diffs: &[E::ScalarField],
    witn_evals: &mut Vec<E::ScalarField>,
) {
    // witn_evals[i] = (evals[i] - fz) / (domain[i]-z)
    for i in 0..domain.size() {
        let num = evals[i] - fz;
        witn_evals.push(num * inv_diffs[i]);
    }
}

/// Evaluate the polynomial given by the evaluations evals over domain at z
/// Assumes that inv_diffs[i] = 1/(domain[i]-z)
/// Note: This assumes that z is not in the evaluation domain
#[inline]
pub fn evaluate_outside<E: Pairing, D: EvaluationDomain<E::ScalarField>>(
    domain: &D,
    evals: &[E::ScalarField],
    z: E::ScalarField,
    inv_diffs: &[E::ScalarField],
) -> E::ScalarField {
    // formula taken from https://dankradfeist.de/ethereum/2021/06/18/pcs-multiproofs.html
    // f(z) = {z^d-1}/d * sum_i (f_i * {w^i}/{z-w^i}), where d is the size of the domain
    let nom = domain.vanishing_polynomial().evaluate(&z);
    let factor = nom / domain.size_as_field_element();
    let mut sum = E::ScalarField::zero();
    for i in 0..domain.size() {
        let term = -domain.element(i) * inv_diffs[i];
        sum += evals[i] * term;
    }
    factor * sum
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_ec::CurveGroup;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{DenseUVPolynomial, Radix2EvaluationDomain};
    use ark_poly::{EvaluationDomain, Evaluations, Polynomial};
    use ark_std::One;
    use ark_std::UniformRand;

    use super::{
        evaluate_outside, find_in_domain, inv_diffs, witness_evals_inside, witness_evals_outside,
    };

    type F = <Bls12_381 as Pairing>::ScalarField;
    type D = Radix2EvaluationDomain<F>;

    /// test function find_in_domain
    #[test]
    fn test_find_in_domain() {
        let mut rng = ark_std::rand::thread_rng();
        let degree = 15;
        let domain = D::new(degree + 1).unwrap();
        let runs = 20;
        for _ in 0..runs {
            // check that it returns None
            // for element outside of domain
            let z = domain.sample_element_outside_domain(&mut rng);
            assert!(find_in_domain::<Bls12_381, D>(&domain, z).is_none());
        }
        // check that it returns the right
        // index for all elements in domain
        for i in 0..domain.size() {
            let fi = find_in_domain::<Bls12_381, D>(&domain, domain.element(i));
            assert!(fi.is_some());
            let j = fi.unwrap();
            assert_eq!(i, j);
        }
    }

    /// test function witness_evals_inside
    #[test]
    fn test_witness_evals_inside() {
        let mut rng = ark_std::rand::thread_rng();
        let degree = 15;
        let domain = D::new(degree + 1).unwrap();
        let runs = 10;
        for _ in 0..runs {
            // sample a random polynomial
            let f = DensePolynomial::rand(domain.size() - 1, &mut rng);
            // evaluate the polynomial on the domain
            let evals = domain.fft(&f.coeffs);
            // test that witness_evals_inside works properly over the entire domain
            for i in 0..domain.size() {
                let z = domain.element(i);
                // compute the witness polynomial by long division
                let fshift = &f - &DensePolynomial::from_coefficients_vec(vec![evals[i]]);
                let div = DensePolynomial::from_coefficients_vec(vec![-z, F::one()]);
                let witness_poly = &fshift / &div;
                let witn_evals_expected = domain.fft(&witness_poly.coeffs);
                // compare with what we get from our function
                let mut witn_evals = Vec::new();
                witness_evals_inside::<Bls12_381, D>(&domain, &evals, i, &mut witn_evals);
                for i in 0..domain.size() {
                    assert_eq!(witn_evals[i], witn_evals_expected[i]);
                }
            }
        }
    }

    /// test function inv_diffs
    #[test]
    fn test_inv_diffs() {
        let mut rng = ark_std::rand::thread_rng();
        let degree = 15;
        let domain = D::new(degree + 1).unwrap();
        let runs = 10;
        for _ in 0..runs {
            // sample an element outside the domain
            let z = domain.sample_element_outside_domain(&mut rng);
            // compute its inverse differences
            let inv_diffs = inv_diffs::<Bls12_381, D>(&domain, z);
            // check that each element is really the inverse
            for i in 0..domain.size() {
                let diff = domain.element(i) - z;
                let prod = diff * inv_diffs[i];
                assert_eq!(prod, F::one());
            }
        }
    }

    /// test function witness_evals_outside
    #[test]
    fn test_witness_evals_outside() {
        let mut rng = ark_std::rand::thread_rng();
        let degree = 15;
        let domain = D::new(degree + 1).unwrap();
        let runs = 10;
        for _ in 0..runs {
            // sample a random polynomial by sampling coefficients
            let f = DensePolynomial::rand(domain.size() - 1, &mut rng);
            // evaluate the polynomial on the domain
            let evals: Vec<F> = domain.fft(&f.coeffs);
            // do a few tests with this polynomial
            for _ in 0..runs {
                // sample a random point outside the domain
                let z = domain.sample_element_outside_domain(&mut rng);
                let inv_diffs = inv_diffs::<Bls12_381, D>(&domain, z);
                let fz = f.evaluate(&z);
                let fshift = &f - &DensePolynomial::from_coefficients_vec(vec![fz]);
                // compute the witness polynomial by long division
                let div = DensePolynomial::from_coefficients_vec(vec![-z, F::one()]);
                let witness_poly = &fshift / &div;
                let witn_evals_expected = domain.fft(&witness_poly.coeffs);
                // compare with what we get from our function
                let mut witn_evals = Vec::new();
                witness_evals_outside::<Bls12_381, D>(
                    &domain,
                    &evals,
                    fz,
                    &inv_diffs,
                    &mut witn_evals,
                );
                for i in 0..domain.size() {
                    assert_eq!(witn_evals[i], witn_evals_expected[i]);
                }
            }
        }
    }

    /// test function evaluate_outside
    #[test]
    fn test_evaluate_outside() {
        let mut rng = ark_std::rand::thread_rng();
        let degree = 15;
        let domain = D::new(degree + 1).unwrap();
        let runs = 20;
        for _ in 0..runs {
            // sample a random polynomial by sampling random evaluations
            let mut evals = Vec::new();
            for _ in 0..domain.size() {
                evals.push(F::rand(&mut rng));
            }
            let evals = Evaluations::from_vec_and_domain(evals, domain);
            let f = evals.interpolate_by_ref();

            // sample a bunch of points and
            // evaluate the polynomial classically
            // and evaluate it using our function
            // result should then be the same
            // tests will most likely be outside of the domain
            for _ in 0..runs {
                let z = F::rand(&mut rng);
                let expected: F = f.evaluate(&z);
                let inv_diffs = inv_diffs::<Bls12_381, D>(&domain, z);
                let obtained: F =
                    evaluate_outside::<Bls12_381, D>(&domain, &evals.evals, z, &inv_diffs);
                assert_eq!(obtained, expected);
            }
        }
    }
}
