package sign

// PARAMETERS
const (
	M               = 8
	N               = 7
	Dbar            = 48
	B               = 430070539612332.205811372782969  // 2^48.61156663661591
	Bsquare         = "184960669042442604975662780477" // B^2
	Kappa           = 23
	LogN            = 8
	SigmaE          = 6.108187070284607
	BoundE          = SigmaE * 2
	SigmaStar       = 172852667880.2713189548230532887787 // 2^37.33075191469097
	BoundStar       = SigmaStar * 2
	SigmaU          = 163961331.5239387
	BoundU          = SigmaU * 2
	KeySize         = 32              // 256 bits
	Q               = 0x1000000004A01 // 48-bit NTT-friendly prime
	QNu             = 0x80000
	QXi             = 0x40000
	TrustedDealerID = 0
	CombinerID      = 1
	Xi              = 30
	Nu              = 29
	EtaEpsilon      = 2.650104
)
