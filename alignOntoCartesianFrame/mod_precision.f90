MODULE mod_precision
  implicit none

  INTEGER(KIND=4),PARAMETER ::           &
          i1 = SELECTED_INT_KIND(1),     &
          i4 = SELECTED_INT_KIND(4),     &
          i8 = SELECTED_INT_KIND(8),     &
          sp = SELECTED_REAL_KIND(6,37), & ! 32-bits precision
          dp = SELECTED_REAL_KIND(15,307)  ! 64-bits precision

  REAL(KIND=8),PARAMETER ::      &
          pi   = 3.14159265359,  &
          g_e  = 2.002319304,    &
          mu_b = 9.27400999E-24, &
          h    = 6.62607015E-34

END MODULE mod_precision
