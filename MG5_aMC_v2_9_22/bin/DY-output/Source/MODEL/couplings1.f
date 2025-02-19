ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_1 = -3.333333D-01*(MDL_EE*MDL_COMPLEXI)
      GC_2 = (2.000000D+00*MDL_EE*MDL_COMPLEXI)/3.000000D+00
      GC_52 = MDL_EE*MDL_COMPLEXI*MDL_QX
      GC_69 = -5.000000D-01*(MDL_CW*MDL_EE*MDL_COMPLEXI)/MDL_SW
      GC_70 = (MDL_CW*MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_SW)
      GC_77 = -1.666667D-01*(MDL_EE*MDL_COMPLEXI*MDL_SW)/MDL_CW
      GC_79 = -((MDL_EE*MDL_COMPLEXI*MDL_QX*MDL_SW)/MDL_CW)
      END
