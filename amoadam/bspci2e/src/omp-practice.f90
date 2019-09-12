program omp_practice
implicit none

integer:: thread_id, nof_threads, omp_get_thread_num, omp_get_num_threads

!$omp parallel private(thread_id)

  thread_id = omp_get_thread_num()
  write(*,*) "This thread is number", thread_id

  if (thread_id .eq. 0) then
    nof_threads = omp_get_num_threads()
    write(*,*) "Number of threads in use is ", nof_threads
  end if

!$omp end parallel

end program omp_practice
