# Written by Henry Gebhardt (2016-2017)

module PkSpectra

export PkSpectrum


using Dierckx  # same as used in the TwoFAST tests


type PkSpectrum
	kk
	pk
	pkspl
	kmin
	kmax
	nslo
	nshi
	kmin_norm
	kmax_norm
end


function PkSpectrum(filename="data/planck_base_plikHM_TTTEEE_lowTEB_lensing_post_BAO_H070p6_JLA_matterpower.dat")
	kk = readdlm(filename)[:,1]
	pk = readdlm(filename)[:,2]
	pkspl = Spline1D(kk, pk)

	# fit low-k using derivative
	k0 = kk[1]
	P0 = pkspl(k0)
	Pp0 = derivative(pkspl, k0)
	nslo = k0 * Pp0 / P0
	kmin_norm = P0 / k0^nslo
	kmin = k0

	# fit high-k using derivative
	k0 = kk[end]
	P0 = pkspl(k0)
	Pp0 = derivative(pkspl, k0)
	nshi = 4 + k0 * Pp0 / P0
	kmax_norm = P0 / (k0^(nshi - 4))
	kmax = k0

	PkSpectrum(kk, pk, pkspl, kmin, kmax, nslo, nshi, kmin_norm, kmax_norm)
end


function pkspectrum(k, pwr)
	if k < pwr.kmin
		return pwr.kmin_norm * k^pwr.nslo
	elseif k > pwr.kmax
		return pwr.kmax_norm * k^(pwr.nshi - 4)
	else
		return pwr.pkspl(k)
	end
end

(pwr::PkSpectrum)(k) = pkspectrum(k, pwr)


end
