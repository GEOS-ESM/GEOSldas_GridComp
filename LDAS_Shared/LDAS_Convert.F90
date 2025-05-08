#include "MAPL_Generic.h"

module LDAS_ConvertMod

  use ESMF
  use MAPL_mod
  use LDAS_DateTimeMod, only: date_time_type
  use LDAS_DateTimeMod, only: get_dofyr_pentad

  implicit none

  private

  public :: esmf2ldas
  public :: string2tile_types

  interface esmf2ldas
     module procedure esmf2ldas_time
  end interface esmf2ldas

contains

  subroutine esmf2ldas_time(esmf_dt, ldas_dt, rc)
    
    type(ESMF_Time),      intent(in)  :: esmf_dt
    type(date_time_type), intent(out) :: ldas_dt
    integer, optional,    intent(out) :: rc

    character(len=*), parameter :: Iam = 'emsf2ldas_time'
    integer :: status

    call ESMF_TimeGet(                                                          &
         esmf_dt,                                                               &
         YY=ldas_dt%year,                                                       &
         MM=ldas_dt%month,                                                      &
         DD=ldas_dt%day,                                                        &
         H=ldas_dt%hour,                                                        &
         M=ldas_dt%min,                                                         &
         S=ldas_dt%sec,                                                         &
         rc=status                                                              &
         )
    VERIFY_(status)

    call get_dofyr_pentad( ldas_dt )

    RETURN_(ESMF_SUCCESS)

  end subroutine esmf2ldas_time

  ! --------------------------------------------
  
  subroutine string2tile_types( string, tile_types)

    ! break string with list of (comma-separated) tile types into vector of strings

    character(len=ESMF_MAXSTR), intent(in)   :: string
    character(10), allocatable, intent(out)  :: tile_types(:)
    
    character(10)  :: outs4(4)
    integer        :: ntype , j, j0

    j  = index(string, ',')      ! identify positions of commas (delimiter)
    ntype  = 1
    j0 = 0

    ! loop through positions of commas
    
    do while (.true.)
       if (j == 0) then
         outs4(ntype) = trim(adjustl(string(j0+1:)))
         exit
       endif
       outs4(ntype) = trim(adjustl(string(j0+1:j0+j-1)))

       j0 = j0+j
       j = index(string(j0+1:), ',')
       ntype = ntype+1
    enddo

    ! assemble output vector of strings
    
    allocate(tile_types(ntype), source=outs4(1:ntype))
    
  end subroutine string2tile_types

  ! --------------------------------------------
  
end module LDAS_ConvertMod

! ========= EOF ===========================================================
