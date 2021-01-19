package bulletproofs

import (
	"crypto/rand"
	"math/big"

	"github.com/ing-bank/zkrp/crypto/p256"
	. "github.com/ing-bank/zkrp/util"
	"github.com/ing-bank/zkrp/util/bn"
)

type MPCContext struct {
	tau1              *big.Int
	tau2              *big.Int
	A                 *p256.P256
	S                 *p256.P256
	T1                *p256.P256
	T2                *p256.P256
	Mu                *big.Int
	Tprime            *big.Int
	InnerProductProof InnerProductProof
	Commit            *p256.P256
	Params            BulletProofSetupParams
}

func StartMPC(params BulletProofSetupParams) (context *MPCContext, tau1H *p256.P256, tau2H *p256.P256) {
	tau1, _ := rand.Int(rand.Reader, ORDER)
	tau2, _ := rand.Int(rand.Reader, ORDER)
	tau1H = new(p256.P256).ScalarMult(params.H, tau1)
	tau2H = new(p256.P256).ScalarMult(params.H, tau2)
	context = &MPCContext{tau1: tau1, tau2: tau2}
	return
}

func (context *MPCContext) PartialProve(secret *big.Int, gamma *big.Int, publicTau1s []*p256.P256, publicTau2s []*p256.P256, alpha *big.Int, rho *big.Int, sL []*big.Int, sR []*big.Int, params BulletProofSetupParams) (taux *big.Int, err error) {
	// aL, aR and commitment: (A, alpha)
	aL, _ := Decompose(secret, 2, params.N) // (41)
	aR, _ := computeAR(aL)                  // (42)
	//alpha, _ := rand.Int(rand.Reader, ORDER)                                   // (43)
	A := commitVector(aL, aR, alpha, params.H, params.Gg, params.Hh, params.N)  // (44)
	S := commitVectorBig(sL, sR, rho, params.H, params.Gg, params.Hh, params.N) // (47)

	// Fiat-Shamir heuristic to compute challenges y and z, corresponds to    (49)
	y, z, _ := HashBP(A, S)

	// compute t1: < aL - z.1^n, y^n . sR > + < sL, y^n . (aR + z . 1^n) >
	vz, _ := VectorCopy(z, params.N)
	vy := powerOf(y, params.N)

	// aL - z.1^n
	naL, _ := VectorConvertToBig(aL, params.N)
	aLmvz, _ := VectorSub(naL, vz)

	// y^n .sR
	ynsR, _ := VectorMul(vy, sR)

	// scalar prod: < aL - z.1^n, y^n . sR >
	sp1, _ := ScalarProduct(aLmvz, ynsR)

	// scalar prod: < sL, y^n . (aR + z . 1^n) >
	naR, _ := VectorConvertToBig(aR, params.N)
	aRzn, _ := VectorAdd(naR, vz)
	ynaRzn, _ := VectorMul(vy, aRzn)

	// Add z^2.2^n to the result
	// z^2 . 2^n
	p2n := powerOf(new(big.Int).SetInt64(2), params.N)
	zsquared := bn.Multiply(z, z)
	z22n, _ := VectorScalarMul(p2n, zsquared)
	ynaRzn, _ = VectorAdd(ynaRzn, z22n)
	sp2, _ := ScalarProduct(sL, ynaRzn)

	// sp1 + sp2
	t1 := big.NewInt(0).Add(sp1, sp2)
	t1 = bn.Mod(t1, ORDER)

	// compute t2: < sL, y^n . sR >
	t2, _ := ScalarProduct(sL, ynsR)
	t2 = bn.Mod(t2, ORDER)

	// compute T1
	T1 := new(p256.P256).ScalarBaseMult(t1) // (53)
	for _, publicTau1 := range publicTau1s {
		T1.Multiply(T1, publicTau1)
	}

	// compute T2
	T2 := new(p256.P256).ScalarBaseMult(t2) // (53)
	for _, publicTau2 := range publicTau2s {
		T2.Multiply(T2, publicTau2)
	}

	// Fiat-Shamir heuristic to compute 'random' challenge x
	x, _, _ := HashBP(T1, T2)

	// compute bl                                                          // (58)
	sLx, _ := VectorScalarMul(sL, x)
	bl, _ := VectorAdd(aLmvz, sLx)

	// compute br                                                          // (59)
	// y^n . ( aR + z.1^n + sR.x )
	sRx, _ := VectorScalarMul(sR, x)
	aRzn, _ = VectorAdd(aRzn, sRx)
	ynaRzn, _ = VectorMul(vy, aRzn)
	// y^n . ( aR + z.1^n sR.x ) + z^2 . 2^n
	br, _ := VectorAdd(ynaRzn, z22n)

	// Compute t` = < bl, br >                                             // (60)
	tprime, _ := ScalarProduct(bl, br)

	// Compute taux = tau2 . x^2 + tau1 . x + z^2 . gamma                  // (61)
	taux = bn.Multiply(context.tau2, bn.Multiply(x, x))
	taux = bn.Add(taux, bn.Multiply(context.tau1, x))
	taux = bn.Add(taux, bn.Multiply(bn.Multiply(z, z), gamma))
	taux = bn.Mod(taux, ORDER)

	// Compute mu = alpha + rho.x                                          // (62)
	mu := bn.Multiply(rho, x)
	mu = bn.Add(mu, alpha)
	mu = bn.Mod(mu, ORDER)

	// Inner Product over (g, h', P.h^-mu, tprime)
	hprime := updateGenerators(params.Hh, y, params.N)

	// SetupInnerProduct Inner Product (Section 4.2)
	var setupErr error
	params.InnerProductParams, setupErr = setupInnerProduct(params.H, params.Gg, hprime, tprime, params.N)
	if setupErr != nil {
		return nil, setupErr
	}
	commit := commitInnerProduct(params.Gg, hprime, bl, br)
	proofip, _ := proveInnerProduct(bl, br, commit, params.InnerProductParams)

	context.A = A
	context.S = S
	context.T1 = T1
	context.T2 = T2
	context.Mu = mu
	context.Tprime = tprime
	context.InnerProductProof = proofip
	context.Commit = commit
	context.Params = params
	return
}

func (context *MPCContext) Aggregate(v *p256.P256, tauxs []*big.Int) (BulletProof, error) {
	// Aggregate taux
	taux := big.NewInt(0)
	for _, partialTaux := range tauxs {
		taux = bn.Add(taux, partialTaux)
	}

	var proof BulletProof
	proof.V = v
	proof.A = context.A
	proof.S = context.S
	proof.T1 = context.T1
	proof.T2 = context.T2
	proof.Taux = taux
	proof.Mu = context.Mu
	proof.Tprime = context.Tprime
	proof.InnerProductProof = context.InnerProductProof
	proof.Commit = context.Commit
	proof.Params = context.Params
	return proof, nil
}
