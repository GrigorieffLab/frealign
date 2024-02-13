      subroutine fftw_plans_2d(n,array,carray,plan_fwd,plan_bwd)
c
      use iso_c_binding
      use fftw33
      implicit none
c
      integer n
      real array(*)
      complex carray(*)
      type(c_ptr) plan_fwd,plan_bwd
c
      plan_fwd=fftwf_plan_dft_r2c_2d(n,n,array,carray,fftw_estimate)
      plan_bwd=fftwf_plan_dft_c2r_2d(n,n,carray,array,fftw_estimate)
c
      return
      end
c
      subroutine fftw_plans_3d(n,array,carray,plan_fwd,plan_bwd)
c
      use iso_c_binding
      use fftw33
      implicit none
c
      integer n
      real array(*)
      complex carray(*)
      type(c_ptr) plan_fwd,plan_bwd
c
      plan_fwd
     . =fftwf_plan_dft_r2c_3d(n,n,n,array,carray,fftw_estimate)
      plan_bwd
     . =fftwf_plan_dft_c2r_3d(n,n,n,carray,array,fftw_estimate)
c
      return
      end
c
      subroutine fftw_fwd(array,carray,plan_fwd)
c
      use iso_c_binding
      use fftw33
      implicit none
c
      real array(*)
      complex carray(*)
      type(c_ptr) plan_fwd,plan_bwd
c
      call fftwf_execute_dft_r2c(plan_fwd,array,carray)
c
      return
      end
c
      subroutine fftw_bwd(array,carray,plan_bwd)
c
      use iso_c_binding
      use fftw33
      implicit none
c
      real array(*)
      complex carray(*)
      type(c_ptr) plan_fwd,plan_bwd
c
      call fftwf_execute_dft_c2r(plan_bwd,carray,array)
c
      return
      end
