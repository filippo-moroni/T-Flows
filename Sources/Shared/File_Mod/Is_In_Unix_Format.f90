!==============================================================================!
  logical function Is_In_Unix_Format(File, file_name)
!------------------------------------------------------------------------------!
!   Find out if file is written in Unix (or Windows) format.                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(File_Type) :: File
  character(len=*) :: file_name
!-----------------------------------[Locals]-----------------------------------!
  integer    :: file_size, file_unit, i
  integer(1) :: byte, next, cnt_13_10
!==============================================================================!

  inquire(file=file_name, size=file_size)

  ! Open the file in binary mode
  call File % Open_For_Reading_Binary(file_name, file_unit, verbose=.false.)

  ! Read a few lines to guess how are new lines being formed
  ! Sequence 13, 10 is used by Windows for line termination
  cnt_13_10 = 0
  do i = 1, min(1024, file_size)  ! don't read beyond the end of file
    read(file_unit) byte
    if(byte .eq. 13) then
      read(file_unit) next
      if(next .eq. 10) then
        cnt_13_10 = cnt_13_10 + 1
      end if
    end if
  end do

  ! Depending on number of occurences of sequence 13, 10, guess the format
  close(file_unit)
  if(cnt_13_10 .lt. 5) then
    Is_In_Unix_Format = .true.
  else
    Is_In_Unix_Format = .false.
  end if

  end function
