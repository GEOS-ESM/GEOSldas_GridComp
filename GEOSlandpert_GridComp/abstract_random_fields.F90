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
