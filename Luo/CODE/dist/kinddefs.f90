module kinddefs

  implicit none
  private

  public :: dp

  integer, parameter :: dp = selected_real_kind(P=15)

end module kinddefs
