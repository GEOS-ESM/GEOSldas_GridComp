module hdf4_fortran_api
  use, intrinsic :: iso_c_binding, only: c_char, c_float, c_int16_t, c_int32_t, c_loc, c_ptr
  use, intrinsic :: iso_fortran_env, only: int16
  implicit none
  private

  public :: hopen, hclose, vfstart, vfend, vsfatch, vsqfnelt, vsfseek, vsfsfld, vsfread, vsfdtch
  public :: sfstart, sfn2index, sfselect, sfginfo, sfrdata, sfendacc, sfend

  interface vsfread
    module procedure vsfread_real
    module procedure vsfread_int16
  end interface

  interface
    function c_hopen(path, path_length, access, ndds) bind(C, name='ldas_h4_hopen') result(status)
      import :: c_char, c_int32_t
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int32_t), value :: path_length, access, ndds
      integer(c_int32_t) :: status
    end function
    function c_vfstart(file_id) bind(C, name='ldas_h4_vfstart') result(status)
      import :: c_int32_t
      integer(c_int32_t), value :: file_id
      integer(c_int32_t) :: status
    end function
    function c_vsfatch(file_id, reference, access, access_length) bind(C, name='ldas_h4_vsfatch') result(status)
      import :: c_int32_t, c_char
      integer(c_int32_t), value :: file_id, reference, access_length
      character(kind=c_char), intent(in) :: access(*)
      integer(c_int32_t) :: status
    end function
    function c_vsqfnelt(vdata_id, records) bind(C, name='ldas_h4_vsqfnelt') result(status)
      import :: c_int32_t
      integer(c_int32_t), value :: vdata_id
      integer(c_int32_t), intent(out) :: records
      integer(c_int32_t) :: status
    end function
    function c_vsfseek(vdata_id, position) bind(C, name='ldas_h4_vsfseek') result(status)
      import :: c_int32_t
      integer(c_int32_t), value :: vdata_id, position
      integer(c_int32_t) :: status
    end function
    function c_vsfsfld(vdata_id, fields, fields_length) bind(C, name='ldas_h4_vsfsfld') result(status)
      import :: c_int32_t, c_char
      integer(c_int32_t), value :: vdata_id, fields_length
      character(kind=c_char), intent(in) :: fields(*)
      integer(c_int32_t) :: status
    end function
    function c_vsfread(vdata_id, data, records, interlace) bind(C, name='ldas_h4_vsfread') result(status)
      import :: c_int32_t, c_ptr
      integer(c_int32_t), value :: vdata_id, records, interlace
      type(c_ptr), value :: data
      integer(c_int32_t) :: status
    end function
    function c_vsfdtch(vdata_id) bind(C, name='ldas_h4_vsfdtch') result(status)
      import :: c_int32_t
      integer(c_int32_t), value :: vdata_id
      integer(c_int32_t) :: status
    end function
    function c_vfend(file_id) bind(C, name='ldas_h4_vfend') result(status)
      import :: c_int32_t
      integer(c_int32_t), value :: file_id
      integer(c_int32_t) :: status
    end function
    function c_hclose(file_id) bind(C, name='ldas_h4_hclose') result(status)
      import :: c_int32_t
      integer(c_int32_t), value :: file_id
      integer(c_int32_t) :: status
    end function
    function c_sfstart(path, path_length, access) bind(C, name='ldas_h4_sfstart') result(status)
      import :: c_char, c_int32_t
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int32_t), value :: path_length, access
      integer(c_int32_t) :: status
    end function
    function c_sfn2index(sd_id, name, name_length) bind(C, name='ldas_h4_sfn2index') result(status)
      import :: c_int32_t, c_char
      integer(c_int32_t), value :: sd_id, name_length
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int32_t) :: status
    end function
    function c_sfselect(sd_id, index) bind(C, name='ldas_h4_sfselect') result(status)
      import :: c_int32_t
      integer(c_int32_t), value :: sd_id, index
      integer(c_int32_t) :: status
    end function
    function c_sfginfo(sds_id, name, name_length, rank, dimsizes, data_type, attributes) bind(C, name='ldas_h4_sfginfo') result(status)
      import :: c_int32_t, c_char
      integer(c_int32_t), value :: sds_id, name_length
      character(kind=c_char), intent(out) :: name(*)
      integer(c_int32_t), intent(out) :: rank, dimsizes(*), data_type, attributes
      integer(c_int32_t) :: status
    end function
    function c_sfrdata(sds_id, start, stride, edge, data) bind(C, name='ldas_h4_sfrdata') result(status)
      import :: c_int32_t, c_ptr
      integer(c_int32_t), value :: sds_id
      integer(c_int32_t), intent(in) :: start(*), stride(*), edge(*)
      type(c_ptr), value :: data
      integer(c_int32_t) :: status
    end function
    function c_sfendacc(sds_id) bind(C, name='ldas_h4_sfendacc') result(status)
      import :: c_int32_t
      integer(c_int32_t), value :: sds_id
      integer(c_int32_t) :: status
    end function
    function c_sfend(sd_id) bind(C, name='ldas_h4_sfend') result(status)
      import :: c_int32_t
      integer(c_int32_t), value :: sd_id
      integer(c_int32_t) :: status
    end function
  end interface

contains
  function hopen(path, access, ndds) result(status)
    character(*), intent(in) :: path
    integer, intent(in) :: access, ndds
    integer :: status
    status = c_hopen(path, len_trim(path), int(access, c_int32_t), int(ndds, c_int32_t))
  end function
  function vfstart(file_id) result(status)
    integer, intent(in) :: file_id
    integer :: status
    status = c_vfstart(int(file_id, c_int32_t))
  end function
  function vsfatch(file_id, reference, access) result(status)
    integer, intent(in) :: file_id, reference
    character(*), intent(in) :: access
    integer :: status
    status = c_vsfatch(int(file_id, c_int32_t), int(reference, c_int32_t), access, len_trim(access))
  end function
  function vsqfnelt(vdata_id, records) result(status)
    integer, intent(in) :: vdata_id
    integer, intent(out) :: records
    integer :: status
    integer(c_int32_t) :: c_records
    status = c_vsqfnelt(int(vdata_id, c_int32_t), c_records)
    records = c_records
  end function
  function vsfseek(vdata_id, position) result(status)
    integer, intent(in) :: vdata_id, position
    integer :: status
    status = c_vsfseek(int(vdata_id, c_int32_t), int(position, c_int32_t))
  end function
  function vsfsfld(vdata_id, fields) result(status)
    integer, intent(in) :: vdata_id
    character(*), intent(in) :: fields
    integer :: status
    status = c_vsfsfld(int(vdata_id, c_int32_t), fields, len_trim(fields))
  end function
  function vsfread_real(vdata_id, data, records, interlace) result(status)
    integer, intent(in) :: vdata_id, records, interlace
    real(c_float), target, intent(inout) :: data(*)
    integer :: status
    status = c_vsfread(int(vdata_id, c_int32_t), c_loc(data(1)), int(records, c_int32_t), int(interlace, c_int32_t))
  end function
  function vsfread_int16(vdata_id, data, records, interlace) result(status)
    integer, intent(in) :: vdata_id, records, interlace
    integer(int16), target, intent(inout) :: data(*)
    integer :: status
    status = c_vsfread(int(vdata_id, c_int32_t), c_loc(data(1)), int(records, c_int32_t), int(interlace, c_int32_t))
  end function
  function vsfdtch(vdata_id) result(status)
    integer, intent(in) :: vdata_id
    integer :: status
    status = c_vsfdtch(int(vdata_id, c_int32_t))
  end function
  function vfend(file_id) result(status)
    integer, intent(in) :: file_id
    integer :: status
    status = c_vfend(int(file_id, c_int32_t))
  end function
  function hclose(file_id) result(status)
    integer, intent(in) :: file_id
    integer :: status
    status = c_hclose(int(file_id, c_int32_t))
  end function
  function sfstart(path, access) result(status)
    character(*), intent(in) :: path
    integer, intent(in) :: access
    integer :: status
    status = c_sfstart(path, len_trim(path), int(access, c_int32_t))
  end function
  function sfn2index(sd_id, name) result(status)
    integer, intent(in) :: sd_id
    character(*), intent(in) :: name
    integer :: status
    status = c_sfn2index(int(sd_id, c_int32_t), name, len_trim(name))
  end function
  function sfselect(sd_id, index) result(status)
    integer, intent(in) :: sd_id, index
    integer :: status
    status = c_sfselect(int(sd_id, c_int32_t), int(index, c_int32_t))
  end function
  function sfginfo(sds_id, name, rank, dimsizes, data_type, attributes) result(status)
    integer, intent(in) :: sds_id
    character(*), intent(out) :: name
    integer, intent(out) :: rank, dimsizes(*), data_type, attributes
    integer :: status
    integer(c_int32_t) :: c_rank, c_dimsizes(2), c_data_type, c_attributes
    name = ' '
    status = c_sfginfo(int(sds_id, c_int32_t), name, len(name), c_rank, c_dimsizes, c_data_type, c_attributes)
    rank = c_rank
    dimsizes(1:2) = c_dimsizes
    data_type = c_data_type
    attributes = c_attributes
  end function
  function sfrdata(sds_id, start, stride, edge, data) result(status)
    integer, intent(in) :: sds_id, start(*), stride(*), edge(*)
    character(kind=c_char), target, intent(inout) :: data(*)
    integer :: status
    integer(c_int32_t) :: c_start(2), c_stride(2), c_edge(2)
    c_start = start(1:2)
    c_stride = stride(1:2)
    c_edge = edge(1:2)
    status = c_sfrdata(int(sds_id, c_int32_t), c_start, c_stride, c_edge, c_loc(data(1)))
  end function
  function sfendacc(sds_id) result(status)
    integer, intent(in) :: sds_id
    integer :: status
    status = c_sfendacc(int(sds_id, c_int32_t))
  end function
  function sfend(sd_id) result(status)
    integer, intent(in) :: sd_id
    integer :: status
    status = c_sfend(int(sd_id, c_int32_t))
  end function
end module hdf4_fortran_api
