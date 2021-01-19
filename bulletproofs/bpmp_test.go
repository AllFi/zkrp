package bulletproofs

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/ing-bank/zkrp/crypto/p256"
	. "github.com/ing-bank/zkrp/util"
	"github.com/ing-bank/zkrp/util/bn"
	"github.com/stretchr/testify/assert"
)

func TestMultiparty(t *testing.T) {
	params, errSetup := Setup(MAX_RANGE_END)
	assert.NoError(t, errSetup)

	value := new(big.Int).SetInt64(int64(300))

	participantsCount := 50
	participantsBlinds := make([]*big.Int, 0)
	for i := 0; i < participantsCount; i++ {
		blind, _ := rand.Int(rand.Reader, ORDER)
		participantsBlinds = append(participantsBlinds, blind)
	}

	participantsContexts := make([]*MPCContext, 0)
	publicTau1s := make([]*p256.P256, 0)
	publicTau2s := make([]*p256.P256, 0)
	for range participantsBlinds {
		context, publicTau1, publicTau2 := StartMPC(params)
		publicTau1s = append(publicTau1s, publicTau1)
		publicTau2s = append(publicTau2s, publicTau2)
		participantsContexts = append(participantsContexts, context)
	}

	alpha, _ := rand.Int(rand.Reader, ORDER)
	rho, _ := rand.Int(rand.Reader, ORDER)
	sL := sampleRandomVector(params.N)
	sR := sampleRandomVector(params.N)

	tauxs := make([]*big.Int, 0)
	for i := 0; i < len(participantsContexts); i++ {
		taux, err := participantsContexts[i].PartialProve(value, participantsBlinds[i], publicTau1s, publicTau2s, alpha, rho, sL, sR, params)
		assert.NoError(t, err)
		tauxs = append(tauxs, taux)
	}

	blindsSum := big.NewInt(0)
	for _, blind := range participantsBlinds {
		blindsSum = bn.Add(blindsSum, blind)
	}
	v, err := CommitG1(value, blindsSum, params.H)
	assert.NoError(t, err)

	proof, err := participantsContexts[0].Aggregate(v, tauxs)
	assert.NoError(t, err)

	ok, err := proof.Verify()
	assert.Equal(t, true, ok)
	assert.NoError(t, err)
}
