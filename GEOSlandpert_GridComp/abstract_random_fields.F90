module abstract_random_fieldsMod
  use nr_ran2_gasdev, only: NRANDSEED
  implicit none
  private
  public :: abstract_random_fields

  type, abstract :: abstract_random_fields
    character(len=:), allocatable :: ID_string
  contains
     procedure(generate_2d_random_field_interface), deferred :: generate_2d_Random_field
     procedure(generate_white_field_interface), deferred :: generate_white_field
     procedure(finalize_random_field_interface), deferred :: finalize
  end type abstract_random_fields

  abstract interface
    subroutine generate_2d_random_field_interface(this, rseed, rfield, rfield2, lx, ly, dx, dy)
      import :: abstract_random_fields, NRANDSEED
      class(abstract_random_fields), intent(inout) :: this
      integer, intent(inout) :: rseed(NRANDSEED)
      real, intent(out) :: rfield(:,:), rfield2(:,:)
      real, intent(in) :: lx, ly, dx, dy
   end subroutine generate_2d_random_field_interface

    subroutine generate_white_field_interface(this, rseed, rfield)
      import :: abstract_random_fields, NRANDSEED
      class(abstract_random_fields), intent(inout) :: this
      integer, intent(inout)  :: rseed(NRANDSEED)
      real, target,intent(out) :: rfield(:,:)
    end subroutine generate_white_field_interface

    subroutine finalize_random_field_interface(this, rc)
      import :: abstract_random_fields
      class(abstract_random_fields), intent(inout) :: this
      integer, optional, intent(out) :: rc
    end subroutine finalize_random_field_interface
  end interface
end module abstract_random_fieldsMod

module StringAbstractRandom_fieldsMapMod
  use abstract_random_fieldsMod

#include "types/key_deferredLengthString.inc"
#define _value class (abstract_random_fields)
#define _value_allocatable
#define _value_equal_defined

! Work around Intel assignment handling for polymorphic map components.
#define _ASSIGN(dest,src) allocate(dest%key,source=src%key); if(allocated(src%value)) allocate(dest%value,source=src%value)
#define _MOVE(dest,src) call move_alloc(from=src%key,to=dest%key); if (allocated(src%value)) call move_alloc(from=src%value,to=dest%value)
#define _FREE(x) deallocate(x%key,x%value)

#define _map StringAbstractRandom_fieldsMap
#define _iterator StringAbstractRandom_fieldsMapIterator
#define _alt

#include "templates/map.inc"

#undef _alt
#undef _iterator
#undef _map
#undef _value
#undef _value_allocatable
#undef _key
#undef _value_equal_defined
#undef _ASSIGN
#undef _MOVE
#undef _FREE
end module StringAbstractRandom_fieldsMapMod
