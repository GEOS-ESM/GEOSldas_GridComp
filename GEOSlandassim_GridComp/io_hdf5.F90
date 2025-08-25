! io_hdf5.F90

! pchakrab, 16 Jan 2014

!---------------------------
! STEPS to read 1D real data
! (For more details, please see the example at the end of this file)
!
! use io_hdf5
! type(hdf5read) :: h5r
! real, allocatable, dimension(:) :: data1D
!
! call h5r%openFile(file)
! for each dset in file
!     call h5r%queryDataset(dsetName, dsetRank, dsetSize)
!     allocate(data1D(dsetSize(1)))
!     call h5r%readDataset(data1D) [READ THE DATASET JUST QUERIED]
! end for
! call h5r%closeFile
! deallocate(data1D)
!---------------------------

module io_hdf5

  use hdf5
  use iso_fortran_env

  implicit none

  private

  integer, parameter :: UNINIT_INT = -99999
  character(len=*), parameter :: UNINIT_STR = ""

  logical, save :: hdf5_inited = .false.
  integer(hid_t), parameter :: INVALID_HID = int(-1,kind=hid_t)

  type, public :: hdf5read
     private
     character(len=1024) :: file_name    = UNINIT_STR
     integer(hid_t)      :: file_id      = INVALID_HID
     character(len=1024) :: dset_name    = UNINIT_STR
     integer(hid_t)      :: dset_id      = INVALID_HID
     integer(hid_t)      :: dspace_id    = INVALID_HID
     integer             :: dset_rank    = 0
     integer(hsize_t)    :: dset_size(7)      = 0_hsize_t
     integer(hsize_t)    :: dset_max_size(7)  = 0_hsize_t
   contains
     ! public
     procedure, public  :: openFile
     procedure, public  :: closeFile
     procedure, public  :: queryDataset
     generic,   public  :: readDataset => readDataset1DReal, readDataset1DReal8, readDataset1DInt, readDataset1DChar24, readDataset2DReal
     ! private
     procedure, private :: readDataset1DReal
     procedure, private :: readDataset1DReal8
     procedure, private :: readDataset1DInt
     procedure, private :: readDataset1DChar24
     procedure, private :: readDataset2DReal
     procedure, private :: uninitDataset
  end type hdf5read

contains

  ! open file
  subroutine openFile(this, filename)
    class (hdf5read), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer :: hdf5err

    this%file_name = filename

    ! Initialize HDF5 fortran interface once per process
    if (.not. hdf5_inited) then
       call h5open_f(hdf5err)
       call checkErrCode_('h5open_f', hdf5err)
       hdf5_inited = .true.
    end if

    call h5fopen_f(this%file_name, H5F_ACC_RDONLY_F, this%file_id, hdf5err)
    call checkErrCode_('h5fopen_f', hdf5err)
  end subroutine openFile

  ! close already opened file
  subroutine closeFile(this)
    class (hdf5read), intent(inout) :: this
    integer :: hdf5err

    if (this%dset_name /= UNINIT_STR) call this%uninitDataset

    if (this%dspace_id >= 0) then
        call h5sclose_f(this%dspace_id, hdf5err)
        if (hdf5err >= 0) this%dspace_id = INVALID_HID
    end if

    if (this%dset_id >= 0) then
        call h5dclose_f(this%dset_id, hdf5err)
        if (hdf5err >= 0) this%dset_id = INVALID_HID
    end if

    if (this%file_id >= 0) then
        call h5fclose_f(this%file_id, hdf5err)
        call checkErrCode_('h5fclose_f', hdf5err)
        this%file_name = UNINIT_STR
        this%file_id   = INVALID_HID
    end if

    ! Do NOT call h5close_f() here.
  end subroutine closeFile

  ! query dataset for number of dims and its shape
  subroutine queryDataset(this, dsetName, dsetRank, dsetSize)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    character(len=*), intent(in) :: dsetName
    integer, intent(out) :: dsetRank
    integer, intent(out) :: dsetSize(7)

    ! local variable
    integer :: hdf5err

    ! ensure that dataset has been uninitialized
    if (this%dset_name/=UNINIT_STR) call this%uninitDataset

    ! set obj param val
    this%dset_name = dsetName

    ! open dataset
    call h5dopen_f(this%file_id, this%dset_name, this%dset_id, hdf5err)
    call checkErrCode_('h5dopen_f', hdf5err)

    ! get data space
    call h5dget_space_f(this%dset_id, this%dspace_id, hdf5err)
    call checkErrCode_('h5dget_space_f', hdf5err)

    ! get number of dims and their sizes
    call h5sget_simple_extent_ndims_f(this%dspace_id, this%dset_rank, hdf5err)
    call checkErrCode_('h5sget_simple_extent_ndims_f', hdf5err)

    call h5sget_simple_extent_dims_f(this%dspace_id, this%dset_size, this%dset_max_size, hdf5err)
    call checkErrCode_('h5sget_simple_extent_dims_f', hdf5err)

    ! return variables
    dsetRank = this%dset_rank
    dsetSize = int(this%dset_size)

  end subroutine queryDataset

  ! uninitialize dataset (close dataset and data space)
  subroutine uninitDataset(this)

    ! input/output variables
    class (hdf5read), intent(inout) :: this

    ! local variable
    integer :: hdf5err

    ! close data space
    if (this%dspace_id >= 0) then
        call h5sclose_f(this%dspace_id, hdf5err)
        call checkErrCode_('h5sclose_f', hdf5err)
        this%dspace_id = INVALID_HID
    end if

    ! close dataset
    if (this%dset_id >= 0) then
        call h5dclose_f(this%dset_id, hdf5err)
        call checkErrCode_('h5dclose_f', hdf5err)
        this%dset_id = INVALID_HID
    end if

    ! uninit obj param val
    this%dset_name = UNINIT_STR
    this%dset_rank = 0
    this%dset_size = 0_hsize_t
    this%dset_max_size = 0_hsize_t

  end subroutine uninitDataset

  ! read 1D character*24 dataset
  subroutine readDataset1DChar24(this, dataChar)
    class (hdf5read), intent(inout) :: this
    character(len=24), intent(out) :: dataChar(:)
    integer :: hdf5err
    integer(hid_t) :: memtype_id

    ! ensure that dset_name is set
    if (this%dset_name == UNINIT_STR) then
       write(*,*) 'ERROR readDataset1DChar24: No open dataset available'
       stop
    end if

    if (this%dset_size(1) == 0) then
       write(*,*) 'Dataset ', trim(this%dset_name), ' is empty'
       return
    end if

    ! Create the memory datatype for the fixed-length strings
    call h5tcopy_f(H5T_FORTRAN_S1, memtype_id, hdf5err)
    call checkErrCode_('h5tcopy_f', hdf5err)

    ! Set size to 24 characters
    call h5tset_size_f(memtype_id, 24_size_t, hdf5err)
    call checkErrCode_('h5tset_size_f', hdf5err)

    ! Read the data using our created type
    call h5dread_f(this%dset_id, memtype_id, dataChar, this%dset_size(1:this%dset_rank), hdf5err)
    call checkErrCode_('h5dread_f', hdf5err)

    ! Close the memory datatype
    call h5tclose_f(memtype_id, hdf5err)
    call checkErrCode_('h5tclose_f', hdf5err)

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DChar24

  ! read 1D real dataset
  subroutine readDataset1DReal(this, data1D)
    class (hdf5read), intent(inout) :: this
    real(REAL32), intent(out) :: data1D(:)
    integer :: hdf5err

    if (this%dset_name==UNINIT_STR) then
       write(*,*) 'ERROR: readDataset1DReal No open dataset available'
       stop
    end if

    if (this%dset_size(1)==0) then
       print *, 'Dataset ', trim(this%dset_name), ' is empty'
    else
       call h5dread_f(this%dset_id, H5T_NATIVE_REAL, data1D, this%dset_size(1:this%dset_rank), hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DReal

  ! read 1D real8 dataset
  subroutine readDataset1DReal8(this, data1D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    real(REAL64), intent(out) :: data1D(:)

    ! local variable
    integer :: hdf5err

    ! check dataset state
    if (this%dset_name == UNINIT_STR) then
       write(*,*) 'ERROR readDataset1DReal8: No open dataset available'
       stop
    end if

    if (this%dset_size(1) == 0) then
       write(*,*) 'Dataset ', trim(this%dset_name), ' is empty'
    else
       call h5dread_f(this%dset_id, H5T_NATIVE_DOUBLE, data1D, this%dset_size(1:this%dset_rank), hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DReal8

  ! read 1D integer dataset
  subroutine readDataset1DInt(this, data1D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    integer, intent(out) :: data1D(:)

    ! local variable
    integer :: hdf5err

    ! check dataset state
    if (this%dset_name == UNINIT_STR) then
        write(*,*) 'ERROR readDataset1DInt: No open dataset available'
        stop
    end if

    if (this%dset_size(1) == 0) then
        write(*,*) 'Dataset ', trim(this%dset_name), ' is empty'
    else
       call h5dread_f(this%dset_id, H5T_NATIVE_INTEGER, data1D, this%dset_size(1:this%dset_rank), hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset1DInt

  ! read 2D real dataset
  subroutine readDataset2DReal(this, data2D)

    ! input/output variables
    class (hdf5read), intent(inout) :: this
    real, intent(out) :: data2D(:,:)

    ! local variable
    integer :: hdf5err

    ! check dataset state
    if (this%dset_name == UNINIT_STR) then
       write(*,*) 'ERROR readDataset2DReal: No open dataset available'
       stop
    end if

    if (this%dset_size(1) == 0) then
       write(*,*) 'Dataset ', trim(this%dset_name), ' is empty'
    else
       call h5dread_f(this%dset_id, H5T_NATIVE_REAL, data2D, this%dset_size(1:this%dset_rank), hdf5err)
       call checkErrCode_('h5dread_f', hdf5err)
    end if

    ! un-initialize dataset just queried/read
    call this%uninitDataset

  end subroutine readDataset2DReal

  ! check return code
  ! (not part of class hdf5read)
  subroutine checkErrCode_(routineName, hdf5errCode)

    ! input/output variables
    character(len=*), intent(in) :: routineName
    integer, intent(in) :: hdf5errCode

    if (hdf5errCode<0) then
       write(*,*) 'ERROR: ', routineName, ' returned NEGATIVE err code: ', hdf5errCode, '. Stopping!'
       stop
    end if

  end subroutine checkErrCode_

end module io_hdf5


! *****************************************************************

#ifdef TEST_IOHDF5

program test_read

  use io_hdf5
  use iso_fortran_env

  implicit none

  character(len=*), parameter :: file_name = '/discover/nobackup/mathomp4/LDAS_Restarts/NGHTLY_TST_TV4000/obs/SMAP/L1C_TB//Y2017/M10/D15/SMAP_L1C_TB_14443_A_20171015T021929_T15160_001.h5'
  character(len=300) :: dsetName

  type(hdf5read) :: h5r
  integer :: dsetRank, dsetSize(7), i

  type myDataType
     real, pointer, dimension(:) :: tb_h_aft => null()     ! aft
     real, pointer, dimension(:) :: lon => null()
     integer, pointer, dimension(:) :: row => null()
     integer, pointer, dimension(:) :: flag => null()
     real(REAL64), pointer, dimension(:) :: tb_time => null()
     character(len=24), pointer, dimension(:) :: tb_time_utc_aft => null()
  end type MyDataType
  type(MyDataType), dimension(1) :: data

  print *, 'HDF5 file: ', trim(file_name)
  print *, ''

  ! open file
  call h5r%openFile(file_name)


  ! query dataset + allocate space + read data
  dsetName = '/Global_Projection/cell_tb_h_aft'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%tb_h_aft(dsetSize(1)))
  call h5r%readDataset(data(1)%tb_h_aft)
  print *, trim(dsetname),'(1:10)'
  print *, data(1)%tb_h_aft(1:10)
  print *, ''

  ! query dataset + allocate space + read data
  dsetName = '/Global_Projection/cell_lon'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%lon(dsetSize(1)))
  call h5r%readDataset(data(1)%lon)
  print *, trim(dsetname),'(1:10)'
  print *, data(1)%lon(1:10)
  print *, ''

  ! query dataset + allocate space + read data
  dsetName = '/Global_Projection/cell_row'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%row(dsetSize(1)))
  call h5r%readDataset(data(1)%row)
  print *, trim(dsetname),'(201:210)'
  print *, data(1)%row(201:210)
  print *, ''

  ! query dataset + allocate space + read data
  dsetName = '/Global_Projection/cell_tb_time_seconds_aft'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%tb_time(dsetSize(1)))
  call h5r%readDataset(data(1)%tb_time)
  print *, trim(dsetname),'(1:10)'
  print *, data(1)%tb_time(1:10)
  print *, ''

  dsetName = '/Global_Projection/cell_tb_time_utc_aft'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%tb_time_utc_aft(dsetSize(1)))
  call h5r%readDataset(data(1)%tb_time_utc_aft)
  print *, trim(dsetname),'(241:250)'
  do i=241,250
     print *, data(1)%tb_time_utc_aft(i)
  end do

  dsetName = '/Global_Projection/cell_tb_qual_flag_h_aft'
  call h5r%queryDataset(dsetName, dsetRank, dsetSize)
  allocate(data(1)%flag(dsetSize(1)))
  call h5r%readDataset(data(1)%flag)
  print *, trim(dsetname),'(241:250)'
  do i=241,250
     print *, data(1)%flag(i)
  end do


  ! close file
  call h5r%closeFile


  ! deallocate memory
  if (associated(data(1)%tb_h_aft)) deallocate(data(1)%tb_h_aft)
  if (associated(data(1)%lon)) deallocate(data(1)%lon)
  if (associated(data(1)%row)) deallocate(data(1)%row)
  if (associated(data(1)%flag)) deallocate(data(1)%flag)
  if (associated(data(1)%tb_time)) deallocate(data(1)%tb_time)
  if (associated(data(1)%tb_time_utc_aft)) deallocate(data(1)%tb_time_utc_aft)

end program test_read

#endif

! =================== EOF ==========================================
