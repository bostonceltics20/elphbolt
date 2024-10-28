program check_fermi_window_adjustment

  use iso_fortran_env, only : r64 => real64, i64 => int64
  use testify_m, only : testify
  use misc, only: print_message, subtitle, timer, exit_with_message
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use electron_module, only: electron
  use phonon_module, only: phonon
  use wannier_module, only: wannier
  use bte_module, only: bte
  use bz_sums, only: calculate_dos, calculate_qTF, calculate_el_dos_fermi, calculate_el_Ws
  use interactions, only: calculate_gReq, calculate_gkRp, calculate_3ph_interaction, &
       calculate_eph_interaction_ibzq, calculate_eph_interaction_ibzk, &
       calculate_echimp_interaction_ibzk, calculate_bound_scatt_rates
      
  implicit none

  integer :: itest
  integer, parameter :: num_tests = 4
  type(testify) :: test_array(num_tests), tests_all

  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(wannier) :: wann
  type(electron) :: el
  type(phonon) :: ph
  type(bte) :: bt
  type(timer) :: t_all, t_event

  !character(:), allocatable :: datalink, curl_arg, workdir, datadir, inputdir

  if(this_image() == 1) then
     write(*, '(A)')  'Regression test on 3C-SiC'
     write(*, '(A, I5)') 'Number of coarray images = ', num_images()
  end if
     
  !Test counter
  itest = 0

  call t_all%start_timer('elphbolt: Fermi-window adjustment')

  call t_event%start_timer('Initialization')

  !Set up crystal
  call crys%initialize

  !Set up numerics data
  call num%initialize(crys)

  !Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, num%qmesh)

  !Test symmetries
  if(this_image() == 1) then
     itest = itest + 1
     test_array(itest) = testify("number of symmetries")
     call test_array(itest)%assert(sym%nsymm, 24_i64)

     itest = itest + 1
     test_array(itest) = testify("symmetry group")
     call test_array(itest)%assert(sym%international, "F-43m")
  end if
  sync all
  !!

  if(num%onlyebte .or. num%drag .or. num%phe &
       .or. num%plot_along_path .or. num%runlevel == 3) then
     !Read EPW Wannier data
     call wann%read(num)

     !Calculate electrons
     call el%initialize(wann, crys, sym, num)

     !Test electron energies
     if(this_image() == 1) then
        itest = itest + 1
        test_array(itest) = testify("electron energies")
        call test_array(itest)%assert( &
             [transpose(el%ens_irred)], &
             [0.1091276329E+02_r64, 0.1459479591E+02_r64, 0.2329596906E+02_r64, 0.2329596906E+02_r64, &
              0.1096209018E+02_r64, 0.1442980731E+02_r64, 0.2334947484E+02_r64, 0.2340698517E+02_r64, &
              0.1079030242E+02_r64, 0.1369997715E+02_r64, 0.2354625098E+02_r64, 0.2354625098E+02_r64, &
              0.1107938027E+02_r64, 0.1467078718E+02_r64, 0.2311065489E+02_r64, 0.2358101936E+02_r64], &
             tol = 1.0e-8_r64)
     end if
     !!
  end if

  call t_event%end_timer('Initialization')

  call t_event%start_timer('Density of states and one-particle scattering rates')

  call subtitle("Calculating density of states...")
  if(num%onlyebte .or. num%drag) then
     !Calculate electron density of states
     call calculate_dos(el, num%tetrahedra)


     !Test electron density of states
     if(this_image() == 1) then
        itest = itest + 1
        test_array(itest) = testify("electron DOS")
        call test_array(itest)%assert( &
             [transpose(el%dos)], &
             [0.6944753689E-02_r64, 0.7200285677E-02_r64, 0.2139544764E-01_r64, 0.2139544797E-01_r64, &
              0.2264851372E-01_r64, 0.3001614228E-02_r64, 0.1790594732E-01_r64, 0.2216848884E-01_r64, &
              0.0000000000E+00_r64, 0.0000000000E+00_r64, 0.6902017130E-02_r64, 0.6902016684E-02_r64, &
              0.9485621719E-02_r64, 0.2888870927E-02_r64, 0.9858175722E-02_r64, 0.1918518823E-02_r64], &
             tol = 1.0e-11_r64)
     end if
    end if
     !!

  if(this_image() == 1) then
     tests_all = testify(test_array)
     call tests_all%report
     
     if(tests_all%get_status() .eqv. .false.) error stop -1
  end if
  sync all
  call t_all%end_timer('elphbolt: Fermi-window adjustment')
end program check_fermi_window_adjustment
