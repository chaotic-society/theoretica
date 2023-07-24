
///
/// @file gauss.h Roots and weights for Gaussian quadrature
///

#ifndef THEORETICA_GAUSS_H
#define THEORETICA_GAUSS_H

namespace theoretica {

	/// @namespace theoretica::tables Tabulated values
	namespace tables {


		// Legendre

		static real legendre_roots_2[2] = {
			-0.577350269189625764506505,
			0.577350269189625764506505
		};

		static real legendre_weights_2[2] = {
			1, 1
		};

		static real legendre_roots_4[4] = {
			-0.861136311594052575194101,
			-0.339981043584856264846549,
			0.339981043584856264738129,
			0.861136311594052575194101
		};

		static real legendre_weights_4[4] = {
			0.347854845137453857491878,
			0.652145154862546142535227,
			0.652145154862546142643647,
			0.347854845137453857491878
		};

		static real legendre_roots_8[8] = {
			-0.960289856497536231715127,
			-0.796666477413626739621007,
			-0.525532409916328985884527,
			-0.183434642495649804990214,
			0.183434642495649804936004,
			0.525532409916328985884527,
			0.796666477413626739621007,
			0.960289856497536231715127
		};

		static real legendre_weights_8[8] = {
			0.101228536290376258828787,
			0.222381034453374470260959,
			0.313706645877887287419484,
			0.362683783378361982948669,
			0.362683783378361982948669,
			0.313706645877887287419484,
			0.222381034453374470260959,
			0.101228536290376258828787
		};

		static real legendre_roots_16[16] = {
			-0.989400934991649929186111,
			-0.944575023073232580589787,
			-0.865631202387831746901999,
			-0.755404408355003026897482,
			-0.617876244402643747097217,
			-0.458016777657227386105678,
			-0.281603550779258913257935,
			-0.0950125098376374401876535,
			0.0950125098376374401605484,
			0.281603550779258913203725,
			0.458016777657227386701989,
			0.617876244402643747097217,
			0.755404408355003026897482,
			0.865631202387831746901999,
			0.944575023073232565953058,
			0.989400934991649929186111
		};

		static real legendre_weights_16[16] = {
			0.0271524594117540951122847,
			0.0622535239386478849447004,
			0.0951585116824928647197876,
			0.124628971255533850669866,
			0.149595988816576744480109,
			0.169156519395002539321907,
			0.182603415044923588682552,
			0.189450610455068496273891,
			0.189450610455068496300996,
			0.182603415044923588967155,
			0.169156519395002538020865,
			0.149595988816576744480109,
			0.124628971255533850669866,
			0.0951585116824928647197876,
			0.0622535239386479161324535,
			0.0271524594117540951122847
		};


		// Laguerre

		static real laguerre_roots_2[2] = {
			0.58578643762690546732,
			3.4142135623730945322
		};

		static real laguerre_weights_2[2] = {
			0.85355339059327000202,
			0.14644660940672634871
		};


		static real laguerre_roots_4[4] = {
			0.32254768961939191742,
			1.7457611011583463777,
			4.5366202969211275316,
			9.3950709123011330332
		};

		static real laguerre_weights_4[4] = {
			0.60315410434164023839,
			0.35741869243780005017,
			0.03888790851500541948,
			0.00053929470556132749637
		};


		static real laguerre_roots_8[8] = {
			0.17027963230510124502,
			0.90370177679937967218,
			2.2510866298661306928,
			4.2667001702876592626,
			7.0459054023934658056,
			10.758516010180995295,
			15.740678641278004193,
			22.863131736889263833
		};

		static real laguerre_weights_8[8] = {
			0.36918858934162849121,
			0.41878678081434484077,
			0.17579498663717180432,
			0.033343492261215595242,
			0.0027945362352256706928,
			9.0765087733582173847e-005,
			8.4857467162725295329e-007,
			1.0480011748715107818e-009
		};


		static real laguerre_roots_16[16] = {
			0.087649410478928189153,
			0.46269632891508138072,
			1.1410577748312272208,
			2.1292836450983810438,
			3.4370866338932063403,
			5.0780186145497663115,
			7.0703385350482545854,
			9.4383143363918924219,
			12.214223368866159374,
			15.441527368781641486,
			19.180156856752924456,
			23.51590569399238859,
			28.578729742881671998,
			34.583398702286930866,
			41.940452647688152288,
			51.701160339543346826
		};

		static real laguerre_weights_16[16] = {
			0.20615171495777392289,
			0.33105785495087120126,
			0.26579577764421131011,
			0.13629693429637680541,
			0.047328928694125499268,
			0.011299900080339313299,
			0.0018490709435266595795,
			0.00020427191530818007586,
			1.4844586873978733152e-005,
			6.8283193308651186878e-007,
			1.8810248410831411265e-008,
			2.8623502429707047297e-010,
			2.1270790332192899185e-012,
			6.2979670025300090354e-015,
			5.050473700032298015e-018,
			4.1614623703730629937e-022
		};

	}

}

#endif