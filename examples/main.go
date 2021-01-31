package main
 
import (
    "fmt"
    "math"
 
    "github.com/ldsec/lattigo/v2/ckks"
)
 
func main() {
 
    params := ckks.DefaultBootstrapSchemeParams[0]
    btpParams := ckks.DefaultBootstrapParams[0]

    fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, logQP = %d, levels = %d, scale= 2^%f, sigma = %f \n", params.LogN(), params.LogSlots(), params.LogQP(), params.Levels(), math.Log2(params.Scale()), params.Sigma())

    kgen := ckks.NewKeyGenerator(params)
    sk, pk := kgen.GenKeyPairSparse(btpParams.H)

    bk := kgen.GenBootstrappingKey(params.LogSlots(), btpParams, sk)

    encoder := ckks.NewEncoder(params)
    decryptor := ckks.NewDecryptor(params, sk)
    encryptor := ckks.NewEncryptorFromPk(params, pk)
    evaluator := ckks.NewEvaluator(params)
    bootstrapper, err := ckks.NewBootstrapper(params, btpParams, bk)
    if err != nil {
        panic(err)
    }
 
   
    data2 := make([]complex128, params.Slots())
    for i := range data2 {
        data2[i] = complex(27.0, 0.0)
    }
    data2[0] = complex(1.0, 1.0)
 
    plaintext2 := encoder.EncodeNew(data2, params.LogSlots())
 
    ciphertext2 := encryptor.EncryptNew(plaintext2)

    printDebug("Initial values", params, ciphertext2, encoder, decryptor)

 
    for ciphertext2.Level() != 10{
        evaluator.DropLevel(ciphertext2, 1)
    }
 
    // bootstrap
    ciphertext2 = bootstrapper.Bootstrapp(ciphertext2)

    printDebug("After Bootstrapp :", params, ciphertext2, encoder, decryptor)

}
 
 
func printDebug(s string, params *ckks.Parameters, ciphertext *ckks.Ciphertext, encoder ckks.Encoder, decryptor ckks.Decryptor){
    values := encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

    fmt.Println()
    fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
    fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
    fmt.Printf("values: %6.10f %6.10f %6.10f %6.10f...\n", values[0], values[1], values[2], values[3])
}