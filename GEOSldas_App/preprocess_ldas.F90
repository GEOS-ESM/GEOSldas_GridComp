
! how to use :
! ./preprocess_ldas option arg1 arg2 arg3

program main

  use preprocess_ldas_routines,     ONLY:    &
       create_mapping,                       &
       createZoominTilefile,                 &
       createZoominBC,                       &
       createZoominVegRestart,               &
       createZoominRestart,                  &
       correctEase,                          &
       convert_pert_rst,                     &
       optimize_latlon
  
  implicit none
  
  character(len=20 ) :: option
  character(len=512) :: arg1
  character(len=512) :: arg2
  character(len=512) :: arg3
  character(len=512) :: arg4
  character(len=512) :: arg5
  character(len=512) :: arg6
  character(len=512) :: arg7
  character(len=512) :: arg8
  character(len=512) :: arg9
  
  character(len=512) :: orig_tile
  character(len=512) :: new_tile
  character(len=512) :: domain_def_file
  character(len=512) :: catch_def_file
  character(len=512) :: out_path 
  character(len=512) :: exp_id 
  character(len=512) :: orig_catch
  character(len=512) :: new_rtm
  character(len=512) :: orig_rtm
  character(len=512) :: new_catch
  character(len=512) :: orig_BC
  character(len=512) :: new_BC
  character(len=512) :: orig_Veg
  character(len=512) :: new_veg
  character(len=512) :: orig_ease
  character(len=512) :: new_ease
  character(len=512) :: f2g_file
  character(len=12 ) :: ymdhm
  character(len=12 ) :: SURFLAY
  character(len=:), allocatable :: new_r, orig_r, tile_types
  integer, allocatable :: int_types(:)
  
  call get_command_argument(1,option)
  call get_command_argument(2,arg1)
  call get_command_argument(3,arg2)
  call get_command_argument(4,arg3)
  call get_command_argument(5,arg4)
  call get_command_argument(6,arg5)
  call get_command_argument(7,arg6)
  call get_command_argument(8,arg7)
  call get_command_argument(9,arg8)
  call get_command_argument(10,arg9)
  
  if( trim(option) == "c_f2g") then

     ! (1) generate 'f2g.txt'
     ! (2) generate tile.domain if it is local

     orig_tile       = arg1
     domain_def_file = arg2
     out_path        = arg3
     catch_def_file  = arg4
     exp_id          = arg5
     ymdhm           = trim(adjustl(arg6))
     SURFLAY         = trim(adjustl(arg7))
     f2g_file        = arg8

     call get_tile_types(trim(arg9), int_types)  

     call create_mapping(orig_tile,domain_def_file,trim(out_path),catch_def_file,trim(exp_id),ymdhm, SURFLAY, f2g_file, int_types)
     
  else if (trim(option) == "zoomin_tile") then
     
     orig_tile  = arg1
     new_tile   = arg2
     f2g_file   = arg3
     call createZoominTilefile(f2g_file, orig_tile,new_tile)
     
  else if (trim(option) == "zoomin_bc" ) then
     
     orig_BC  = arg1
     new_BC   = arg2
     f2g_file = arg3

     call createZoominBC(f2g_file, orig_BC, new_BC)
     
  else if (trim(option) == "zoomin_vegrst") then

     orig_veg = arg1
     new_veg  = arg2
     f2g_file = arg3

     call  createZoominVegRestart(f2g_file, orig_veg, new_veg)      

  else if (trim(option) == "zoomin_mwrtmrst") then

     orig_rtm = arg1
     new_rtm  = arg2
     f2g_file = arg3

     call  createZoominRestart(f2g_file, orig_rtm, new_rtm, 100)      

  else if (trim(option) == "zoomin_catchrst") then

     orig_catch = arg1
     new_catch  = arg2
     f2g_file   = arg3

     call createZoominRestart(f2g_file, orig_catch, new_catch, 100)

  else if (trim(option) == "zoomin_landicerst") then

     orig_r = trim(arg1)
     new_r  = trim(arg2)
     f2g_file   = trim(arg3)

     call createZoominRestart(f2g_file, orig_r, new_r, 20)

  else if (trim(option)=="correctease") then

     orig_ease = arg1 
     new_ease  = arg2

     call correctEase(orig_ease,new_ease) 

  else if (trim(option)=="c_convert_pert") then

     out_path  = arg3
     exp_id    = arg4

     call convert_pert_rst(arg1,arg2, out_path,exp_id)
     
  else if (trim(option) == "optimize") then
     
     call get_tile_types(trim(arg5), int_types) 
     call optimize_latlon(arg1,arg2, arg3, arg4, int_types)
     
  else
     
     print*, " wrong preprocess option:",option
     
  end if

contains

  subroutine get_tile_types(str_types, int_types)
    character(*), intent(in) :: str_types
    integer, allocatable, intent(out) :: int_types(:)
    integer :: n, Length, from, to, i
    n = 1
    Length = len(str_types)
    do i = 1, Length
      if (str_types(i:i) == '_') n = n+1
    enddo
    allocate(int_types(n))
    from = 0
    to   = 1
    n    = 1
    do while (to <= Length)
      if (str_types(to:to) == "_") then
         read (unit=str_types(from+1:to-1),fmt=*) int_types(n)
         n = n + 1
         from = to
      endif
      to = to + 1
    enddo
    read (unit=tile_types(from+1:to-1),fmt=*) int_types(n)
  end subroutine get_tile_types 

end program main

! ====================== EOF =======================================================
