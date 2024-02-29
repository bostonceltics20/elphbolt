module resource_module
  !Module containing the data type for resource management.

#ifdef _OPENACC
  use openacc
#endif
  use params, only: r64, i64
  use misc, only: create_set
  
  implicit none
  
  public resource
  
  type resource
     integer :: num_cpus
     !! Number of lonely cpus (not a gpu manager)
     integer :: num_gpus
     !! Number of gpus
     integer :: this_node
     !! ID of node where this image resides
     logical, allocatable :: gpu_manager[:]
     !! Is this image a gpu manager?
     character*(30) :: gpu_name = trim("None"), &
          gpu_vendor = trim("None"), gpu_driver = trim("None")
     !! Information about gpu device
     
   contains
     procedure :: initialize, report, balance_load
  end type resource

contains
  
  subroutine initialize(self)
    class(resource), intent(out) :: self

    integer :: igpus, im, iset
    integer :: devicenum
#ifdef _OPENACC
    integer(acc_device_property):: property
#endif

    character*(100) :: string
    character(len=10), allocatable :: hostname(:)[:]
    character(len=10), allocatable :: hostname_set(:)

    allocate(self%gpu_manager[*])
        
    allocate(hostname(num_images())[*])
    call hostnm(hostname(this_image()))
    if(this_image() == 1) then
       do im = 1, num_images()
          hostname(im) = hostname(im)[im]
       end do
       !print*, 'hostname: ', hostname(:)
    end if
    sync all
    call co_broadcast(hostname, 1)
    sync all
    
    call create_set(hostname, hostname_set)
    
    self%num_gpus = 0
#ifdef _OPENACC
    !Bug? The following returns 1 from all images, even if I throw multiple
    !gpu-equipped nodes at the program. 
    !self%num_gpus = acc_get_num_devices(acc_device_default)
    !
    !A way around?
    self%num_gpus = size(hostname_set)
    !print*, 'Num devices: ', acc_get_num_devices(acc_device_default)
#endif
    self%num_cpus = num_images() - self%num_gpus

    self%this_node = findloc(hostname_set, hostname(this_image()), 1)
    self%gpu_manager = &
         (this_image() == findloc(hostname, hostname_set(self%this_node), 1)) &
         .and. self%num_gpus > 0
    
#ifdef _OPENACC
    if(self%gpu_manager) then
       igpus = 0
       property = acc_property_name
       call acc_get_property_string(igpus, acc_get_device_type(), &
            property, string)
       self%gpu_name = trim(string)

       property = acc_property_vendor
       call acc_get_property_string(igpus, acc_get_device_type(), &
            property, string)
       self%gpu_vendor = trim(string)

       property = acc_property_driver
       call acc_get_property_string(igpus, acc_get_device_type(), &
            property, string)
       self%gpu_driver = trim(string)
       
!!$       do igpus = 0, self%num_gpus - 1 !Mind the 0 based indexing of openacc
!!$          property = acc_property_name
!!$          call acc_get_property_string(igpus, acc_get_device_type(), &
!!$               property, string)
!!$          self%gpu_name = trim(string)
!!$
!!$          property = acc_property_vendor
!!$          call acc_get_property_string(igpus, acc_get_device_type(), &
!!$               property, string)
!!$          self%gpu_vendor = trim(string)
!!$
!!$          property = acc_property_driver
!!$          call acc_get_property_string(igpus, acc_get_device_type(), &
!!$               property, string)
!!$          self%gpu_driver = trim(string)
!!$       end do
    end if
#endif
  end subroutine initialize

  subroutine balance_load(self, split, total_load, chunk, istart, iend, num_active_images)
    !! A simple gpus/cpus load balancer.
    
    class(resource), intent(in) :: self
    real(r64), intent(in) :: split
    integer(i64), intent(in) :: total_load
    integer(i64), intent(out) :: chunk[*], istart[*], iend[*], num_active_images

    !Locals
    integer(i64) :: im, offset

    integer(i64) :: gpu_load, cpu_load, load_per_cpu, load_per_gpu, remainder
    integer(i64), allocatable :: activate[:]

    allocate(activate[*])
    activate = 0
    sync all

    if(this_image() == 1) then       
       !Divide irred phonon states among images/gpu
       gpu_load = nint(split*total_load)
       cpu_load = total_load - gpu_load

       !Maximum Number of active images
       num_active_images = min(gpu_load + cpu_load, num_images())

       !Approximate number of points per image
!!$       load_per_gpu = gpu_load/self%num_gpus
!!$       load_per_cpu = cpu_load/(num_active_images - self%num_gpus)

       load_per_gpu = ceiling(1.0*gpu_load/self%num_gpus)
       load_per_cpu = ceiling(1.0*cpu_load/(num_active_images - self%num_gpus))

!!$       print*, total_load, self%num_gpus, cpu_load, gpu_load, load_per_cpu, load_per_gpu, &
!!$            modulo(cpu_load, num_active_images - self%num_gpus)

       remainder = total_load - load_per_gpu*self%num_gpus - load_per_cpu*self%num_cpus
       
       !print*, remainder
       
       offset = 0
       do im = 1, num_active_images
          if(self%gpu_manager[im]) then
             chunk[im] = min(load_per_gpu, total_load - offset)
          else
             chunk[im] = min(load_per_cpu, total_load - offset)
          end if

          if(chunk[im] > 0) then
             istart[im] = offset + 1
             iend[im] = istart[im] + chunk[im] - 1
             offset = iend[im]
             activate[im] = 1
          else
             activate[im] = 0
             istart[im] = 0
             iend[im] = 0
          end if
          
!!$          print*, 'image, gpu_manager, chunk, istart, iend:', &
!!$               im, self%gpu_manager[im], chunk[im], istart[im], iend[im]
       end do

!!$       if(remainder > 0) then
!!$          !Put the remainder load in the last image
!!$          chunk[num_active_images] = chunk[num_active_images] + remainder
!!$          iend[num_active_images] = iend[num_active_images] + remainder
!!$
!!$          print*, 'Updated! image, gpu_manager, chunk, istart, iend:', &
!!$               num_active_images, self%gpu_manager[num_active_images], &
!!$               chunk[num_active_images], istart[num_active_images], &
!!$               iend[num_active_images]
!!$       end if
    end if
    sync all
    call co_sum(activate)
    sync all

    !Update number of active images
    num_active_images = activate
  end subroutine balance_load
  
  subroutine report(self)
    class(resource), intent(in) :: self

    if(self%gpu_manager) then
       print*, "Node, image, device name, device vendor, device driver", &
            self%this_node, this_image(), self%gpu_name, &
            self%gpu_vendor, self%gpu_driver
    end if
  end subroutine report
end module resource_module
