	program main
	  call imax_verify_(72)
	  call imax_verify_(144)
	  call imax_verify_(288)
	  call imax_verify_(360)
	  call imax_verify_(540)
	  call imax_verify_(576)
	contains
	include 'imax_verify.h'
	end program main
