#include <complex>
#include <cmath>

namespace fft
{
	typedef long long LONG;
	//constexpr double M_PI = 3.1415926535897932384626433;
		
	template <typename Tn>
	std::complex<Tn>* fft(const std::complex<Tn>* data, std::complex<Tn>* F,
		LONG len, LONG start = -1, LONG end = -1, LONG step = -1)
	{
		if (!data || !F)
		{
			return NULL;
		}

		if (start == -1)
		{
			start = 0;
			end = len;
			step = 1;
		}

		if (len == 1)
		{
			F[start] = data[start];
			return F;
		}
		if ((len & (len-1)))
			// len is not an interger power of 2
		{
			return NULL;
		}

		fft(data, F, len / 2, start, end, step * 2);
		// calculate the even part.
		fft(data, F, len / 2, start + step, end, step * 2);
		// calculate the odd part.


		double theta = -2 * M_PI / len;
		double theta0 = -2 * M_PI / len;
		F[start] += F[start + step];
		std::complex<Tn> omega(std::cos(theta), std::sin(theta));
		auto omega0 = omega;
		for (LONG i = 1, j = start+step, k = start+2*step; i < len / 2; i += 1, j+=step, k+=2*step)
		{
			//std::complex<Tn> omega(std::cos(theta), std::sin(theta));
			F[j] = F[k]		// the even part, even[i]
				+ omega * F[k+step];	// the odd part, odd[i]
			omega *= omega0;
		}

		//std::complex<Tn> omega0(std::cos(theta0), std::sin(theta0));
		//auto omega = omega0;
		//for (LONG i = start+step, j = start+2*step, k = 1; k<len/2; i+=step, j+=2*step, k+=1, theta+=theta0)
		//{
		//	//std::complex<Tn> omega(std::cos(theta), std::sin(theta));
		//	F[i] += F[j];
		//	F[i] += omega * F[j + step];
		//	omega *= omega0;
		//}

		for (LONG i = start+step, j = start+(len-1)*step, k = 1; k < len / 2;
			i += step, j-=step, k+=1)
		{
			F[j] = std::complex<Tn>(F[i].real(), -F[i].imag());

		}
		{
			const auto N2 = start + (len / 2) * step;
			F[N2] = 0;
			for (LONG i = start; i < end; i += step * 2)
			{
				F[N2] += data[i];
				F[N2] -= data[i + step];
			}
		}
		return F;
	}


};

extern "C"
_declspec(dllexport)
double* fft_C(double* data, double* result, __int64 len)
// This function is a C-like style of fft::fft() function for external usage.
// data must have 2*len length. Its format is [real 0, imag 0, real 1, imag 1, ... ]. result is the same.
{
	if ((len & (len - 1)))
	{
		return NULL;
	}
	std::complex<double>* ddata = reinterpret_cast<std::complex<double>*> (data);
	std::complex<double>* rresult = reinterpret_cast<std::complex<double>*> (result);


	return reinterpret_cast<double*>(fft::fft(ddata, rresult, len));
}