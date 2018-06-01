# This module, "PerformanceStats", provides the necessary machinery to use
# @timed efficiently. That is, with this module you can use the "+=" operator
# on the object returned by the "@timed" macro. That is, you can write the
# following code:
# 
# 	using PerformanceStats
#	stats = @timed res = somefunction()
#	stats += @timed red = someotherfunction()
#	timed_println("resources: ", stats)

module PerformanceStats

export +, timed_println
import Base.+  # This is the function we want to extend

if VERSION <= v"0.6.9"
	using Compat
	import Base.round
	round(n; digits=0) = round(n, digits)
end


function +(x::Base.GC_Diff, y::Base.GC_Diff)
	return Base.GC_Diff(
		x.allocd	+ y.allocd,	# Bytes allocated
		x.malloc	+ y.malloc,	# Number of GC aware malloc()
		x.realloc	+ y.realloc,	# Number of GC aware realloc()
		x.poolalloc	+ y.poolalloc,	# Number of pool allocations
		x.bigalloc	+ y.bigalloc,	# Number of big (non-pool) allocations
		x.freecall	+ y.freecall,	# Number of GC aware free()
		x.total_time	+ y.total_time,	# Time spent in garbage collection
		x.pause		+ y.pause,	# Number of GC pauses
		x.full_sweep	+ y.full_sweep)	# Number of GC full collection
end

function +(x::Tuple{Any, Float64, Int64, Float64, Base.GC_Diff},
	   y::Tuple{Any, Float64, Int64, Float64, Base.GC_Diff})
	return (Nothing,	# value of the expression
		x[2] + y[2],	# elapsed time
		x[3] + y[3],	# total bytes allocated
		x[4] + y[4],	# garbage collection time
		x[5] + y[5])	# GC_Diff object
end

function prettyprint_getunits_now(value, units, factor)
	if value == 0 || value == 1
		return value, units[1]
	end
	u = ceil(Int, log(value) / log(factor))
	u = min(length(units), u)
	if u == 1
		return value, units[u]
	end
	number = value / factor^(u - 1)
	return round(number, digits=3), units[u]
end

function timed_println(prefix, x::Tuple{Any, Float64, Int64, Float64, Base.GC_Diff})
	mem_units = [" byte", " KB", " MB", " GB", " TB", " PB"]
	cnt_units = ["", " k", " M", " G", " T", " P"]
	allocs = x[5].malloc + x[5].realloc + x[5].poolalloc + x[5].bigalloc
	bytes, mb = prettyprint_getunits_now(x[3], mem_units, 1024)
	allocs, ma = prettyprint_getunits_now(allocs, cnt_units, 1000)
	s = string("$(round(x[2],digits=6)) sec ($(allocs)$(ma) allocations: $(bytes)$(mb), $(round(100*x[4]/x[2],digits=2))% gc time)")
	println(prefix, s)
end

end
