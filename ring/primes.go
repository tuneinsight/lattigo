package ring

// Pi60 are the first hundred (from 0x800000000000000 and upward) 60bit Primes allowing NTT for N = 65536
var Pi60 = []uint64{576460752308273153, 576460752315482113, 576460752319021057, 576460752319414273, 576460752321642497,
	576460752325705729, 576460752328327169, 576460752329113601, 576460752329506817, 576460752329900033,
	576460752331210753, 576460752337502209, 576460752340123649, 576460752342876161, 576460752347201537,
	576460752347332609, 576460752352837633, 576460752354017281, 576460752355065857, 576460752355459073,
	576460752358604801, 576460752364240897, 576460752368435201, 576460752371187713, 576460752373547009,
	576460752374333441, 576460752376692737, 576460752378003457, 576460752378396673, 576460752380755969,
	576460752381411329, 576460752386129921, 576460752395173889, 576460752395960321, 576460752396091393,
	576460752396484609, 576460752399106049, 576460752405135361, 576460752405921793, 576460752409722881,
	576460752410116097, 576460752411033601, 576460752412082177, 576460752416145409, 576460752416931841,
	576460752421257217, 576460752427548673, 576460752429514753, 576460752435281921, 576460752437248001,
	576460752438558721, 576460752441966593, 576460752449044481, 576460752451141633, 576460752451534849,
	576460752462938113, 576460752465952769, 576460752468705281, 576460752469491713, 576460752472375297,
	576460752473948161, 576460752475389953, 576460752480894977, 576460752483254273, 576460752484827137,
	576460752486793217, 576460752486924289, 576460752492691457, 576460752498589697, 576460752498720769,
	576460752499507201, 576460752504225793, 576460752505405441, 576460752507240449, 576460752507764737,
	576460752509206529, 576460752510124033, 576460752510779393, 576460752511959041, 576460752514449409,
	576460752516284417, 576460752519168001, 576460752520347649, 576460752520609793, 576460752522969089,
	576460752523100161, 576460752524279809, 576460752525852673, 576460752526245889, 576460752526508033,
	576460752532013057, 576460752545120257, 576460752550100993, 576460752551804929, 576460752567402497,
	576460752568975361, 576460752573431809, 576460752580902913, 576460752585490433, 576460752586407937}

// Qi60 are the last hundred (from 0xfffffffffffffff and downward) 60bit Primes allowing NTT for N = 65536
var Qi60 = []uint64{1152921504606584833, 1152921504598720513, 1152921504592429057, 1152921504581419009, 1152921504580894721,
	1152921504578273281, 1152921504577748993, 1152921504577486849, 1152921504568836097, 1152921504565166081,
	1152921504563331073, 1152921504556515329, 1152921504555466753, 1152921504554156033, 1152921504552583169,
	1152921504542883841, 1152921504538951681, 1152921504537378817, 1152921504531873793, 1152921504521650177,
	1152921504509853697, 1152921504508280833, 1152921504506970113, 1152921504495697921, 1152921504491241473,
	1152921504488620033, 1152921504479444993, 1152921504470794241, 1152921504468172801, 1152921504462929921,
	1152921504462667777, 1152921504455589889, 1152921504447987713, 1152921504442482689, 1152921504436191233,
	1152921504427278337, 1152921504419414017, 1152921504409190401, 1152921504403947521, 1152921504396869633,
	1152921504395821057, 1152921504373014529, 1152921504369344513, 1152921504368558081, 1152921504364625921,
	1152921504362790913, 1152921504361218049, 1152921504353615873, 1152921504337887233, 1152921504337625089,
	1152921504321372161, 1152921504314032129, 1152921504303022081, 1152921504301449217, 1152921504288342017,
	1152921504287293441, 1152921504286769153, 1152921504282836993, 1152921504274972673, 1152921504266321921,
	1152921504256622593, 1152921504253739009, 1152921504245088257, 1152921504241942529, 1152921504240107521,
	1152921504239583233, 1152921504238010369, 1152921504234078209, 1152921504231718913, 1152921504230670337,
	1152921504227524609, 1152921504214417409, 1152921504207339521, 1152921504205504513, 1152921504204193793,
	1152921504190824449, 1152921504179552257, 1152921504177192961, 1152921504176668673, 1152921504174309377,
	1152921504172474369, 1152921504164872193, 1152921504162512897, 1152921504139706369, 1152921504134987777,
	1152921504132628481, 1152921504122142721, 1152921504120832001, 1152921504116899841, 1152921504105627649,
	1152921504101957633, 1152921504100384769, 1152921504096452609, 1152921504093306881, 1152921504078364673,
	1152921504067092481, 1152921504066306049, 1152921504057917441, 1152921504053723137, 1152921504050839553}
