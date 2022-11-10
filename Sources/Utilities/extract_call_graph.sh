#!/bin/bash

BLACK='\U001B[30m'
RED='\U001B[31m'
GREEN='\U001B[32m'
YELLOW='\U001B[33m'
BLUE='\U001B[34m'
MAGENTA='\U001B[35m'
CYAN='\U001B[36m'
WHITE='\U001B[37m'
RESET='\U001B[0m'

LIGHT_BLACK='\U001B[30;1m'
LIGHT_RED='\U001B[31;1m'
LIGHT_GREEN='\U001B[32;1m'
LIGHT_YELLOW='\U001B[33;1m'
LIGHT_BLUE='\U001B[34;1m'
LIGHT_MAGENTA='\U001B[35;1m'
LIGHT_CYAN='\U001B[36;1m'
LIGHT_WHITE='\U001B[37;1m'
RESET='\U001B[0m'

tabs 60

#------------------------------------------------------------------------------#
# Print_usage
#------------------------------------------------------------------------------#
print_usage() {
  echo "#======================================================================"
  echo "# Utility for extract call graphs from T-Flows"
  echo "#----------------------------------------------------------------------"
  echo "# Proper usage: "
  echo "#"
  echo "# ./Utilities/extract_call_graph.sh <Source> [-f <Force_Dir>]"
  echo "#"
  echo "# where Source  is the procedure name for which you want to perform"
  echo "# the analysis, such as: Main_Con, Compute_Energy, Main_Div, hence"
  echo "# case sensitive, without the .f90 extension."
  echo "#"
  echo "# In cases where the same module name is used in more than one direc-"
  echo "# tory, you can use the second (always -e) and the third argument to "
  echo "# exclude some directories from the search.  At the time of writing "
  echo "# this, Allocate_Cells is defined in Grid_Mod and Refines_Mod, but ."
  echo "# there could be more such examples."
  echo "#"
  echo "# NOTE: The script is supposed to be executed from: T-Flows/Sources!"
  echo "#----------------------------------------------------------------------"
}

glo_procedure=""
glo_module=""

#------------------------------------------------------------------------------#
# Sets $glo_procedure and $glo_module
#------------------------------------------------------------------------------#
extract_procedure_and_module() {

  full_name=$1

  # It rarely gets cooler than this: extract the calling function pattern :-)
  pattern="${full_name//[^%^(^)]}"

  # This should take care of Comm_Mod_Friendly_Function
  if [[ "$pattern" == "" || "$full_name" == *"_Mod_"* ]]; then
    glo_module=$(echo $full_name | awk -F '_Mod_' '{print $1}')"_Mod"
    glo_procedure=$(echo $full_name | awk -F '_Mod_' '{print $2}')
    glo_procedure=$(cut -d \( -f 1 <<< $glo_procedure)

  # Take care of cases such as Grid%Comm%Sendrecv_Real_Arrays
  elif [[ "$pattern" == "%%" ]]; then
    glo_module=$(cut -d % -f 2 <<< $full_name)
    glo_module=$(cut -d \( -f 1 <<< $glo_module)
    glo_procedure=$(cut -d % -f 3 <<< $full_name)
    glo_procedure=$(cut -d \( -f 1 <<< $glo_procedure)

  # Likes of Profiler%Start() and Grid(d)%Calculate() ...
  # Profiler%Update_By_Rank(Profiler%previously_running)
  elif [[ "$pattern" == "%()"       ||     \
          "$pattern" == "%("        ||     \
          "$pattern" == "%(())"     ||     \
          "$pattern" == "()%()"     ||     \
          "$pattern" == "()%("      ||     \
          "$pattern" == "%(%)"      ||     \
          "$pattern" == "%(()"      ||     \
          "$pattern" == "%(%%"      ||     \
          "$pattern" == "%(%%)"     ||     \
          "$pattern" == "%(%())"    ||     \
          "$pattern" == "%(%%())"   ||     \
          "$pattern" == "%((%)"     ||     \
          "$pattern" == "%((%))"    ||     \
          "$pattern" == "%(%%()"    ||     \
          "$pattern" == "%(%%())"   ||     \
          "$pattern" == "%(%()%()"  ||     \
          "$pattern" == "%(()%)" ]]; then
    glo_module=$(cut -d %  -f 1 <<< $full_name)
    glo_module=$(cut -d \( -f 1 <<< $glo_module)
    glo_procedure=$(cut -d %  -f 2 <<< $full_name)
    glo_procedure=$(cut -d \( -f 1 <<< $glo_procedure)

  # This is for global functions and external functions like Probe_1d()
  # or system_clock(count_rate=Profiler%sys_count_rate)
  elif [[ "$pattern" == "(())"      ||       \
          "$pattern" == "()"        ||       \
          "$pattern" == "("         ||       \
          "$pattern" == ""          ||       \
          "$pattern" == "(%()%())"  ||       \
          "$pattern" == "(%()%)"    ||       \
          "$pattern" == "(()())"    ||       \
          "$pattern" == "(%)" ]]; then
    glo_module="" # empty string (maybe none or Shared?)
    glo_procedure=$(cut -d \( -f 1 <<< $full_name)

  # or Comm_Mod_Global_Sum_Real(Profiler%funct_time(i_fun))
  elif [[ "$pattern" == "(%())" ]]; then
    glo_module=$(echo $full_name | awk -F '_Mod_' '{print $1}')"_Mod"
    glo_procedure=$(echo $full_name | awk -F '_Mod_' '{print $2}')
    glo_procedure=$(cut -d \( -f 1 <<< $glo_procedure)

  elif [[ "$pattern" == ")%()" ]]; then
    echo "Warning: complex pattern 'call' is not in the beginning"
    echo "Full name is: " $full_name

  else
    echo "Unknown pattern: ", $pattern
    echo "in             : ", $full_name
    exit
  fi
}

#------------------------------------------------------------------------------#
# Browse through all directories looking for a call graph
#------------------------------------------------------------------------------#
extract_call_graph() {

  #-----------------------
  #   Handle parameters
  #-----------------------

  # First parameter is the procedure name you seek
  procedure_name_you_seek="$1"
  procedure_file_you_seek="$procedure_name_you_seek"".f90"

  # Second parameter is the module in which the procedure resides
  local module_in_which_you_seek="$2"

  local next_level=`expr $next_level + 1`
  local this_level=`expr $next_level - 1`

  #----------------------------------------------#
  #   Get the full path of the module you seek   #
  #----------------------------------------------#
  if [ $module_in_which_you_seek ] && [ $exclude_dir ]; then
    local full_path_you_seek=$(find . -name $procedure_file_you_seek | grep $module_in_which_you_seek | grep -v $exclude_dir)
  elif [ $exclude_dir ]; then
    local full_path_you_seek=$(find . -name $procedure_file_you_seek | grep -v $exclude_dir)
  elif [ $module_in_which_you_seek ]; then
    local full_path_you_seek=$(find . -name $procedure_file_you_seek | grep $module_in_which_you_seek)
  else
    local full_path_you_seek=$(find . -name $procedure_file_you_seek)
  fi

  # This command counts number of occurrences of modules name in the result of
  # command find. If it is more than one, the same file is in more directories
  n=$(echo "$full_path_you_seek" | tr " " "\n" | grep -c "$procedure_name_you_seek")
  if [ $n -gt 1 ]; then
    echo "Ambiguity: procedure "$procedure_name_you_seek" found in more than one directory, here is the list:"
    for path in ${full_path_you_seek[*]}; do
      echo $path
    done
    echo "Exclude all but one directory with the command line argument -e <directory>"
    exit
  fi

  # If the definition of this procedure is found, carry on
  if [ $full_path_you_seek ]; then

    #-----------------------------------------------------
    #   Storing results of the grep command in an array
    #-----------------------------------------------------
    local called_procedures=($(grep '\ \ call\ ' $full_path_you_seek | awk '{print $2$3$4$5$6$7$8$9}' | tr -d ,))

    local called_modules=$called_procedures   # just to declare

    # echo "FETCHED PROCEDURES:"
    # Correct the lines in which call was not the first comman, something
    # like: if(some_condition_is_met) call The_Function_You_Seek
    for proc in "${!called_procedures[@]}"; do
      first_four=${called_procedures[proc]:0:4}
      if [[ "$first_four" == "call" ]]; then   # call was not the first statement in the line
        called_procedures[proc]=$(awk -F'call ?' '{print $2}' <<< ${called_procedures[proc]})
      fi
    done

    # At this point, $called procedures has a form like: Profiler%Start('Main')
    # From this mess, extract the module name and the procedure element wise
    for proc in "${!called_procedures[@]}"; do
      extract_procedure_and_module ${called_procedures[proc]}  # sets $glo_procedure and $glo_module

      # Typical substitues:
      if [ ! "$glo_module" == "" ]; then
        if [   "$glo_module" == "Msg" ];       then glo_module="Message_Mod";   fi
        if [   "$glo_module" == "Message" ];   then glo_module="Message_Mod";   fi
        if [   "$glo_module" == "Tok" ];       then glo_module="Tokenizer_Mod"; fi
        if [   "$glo_module" == "Line" ];      then glo_module="Tokenizer_Mod"; fi
        if [   "$glo_module" == "Vof" ];       then glo_module="Vof_Mod";       fi
        if [   "$glo_module" == "Grid" ];      then glo_module="Grid_Mod";      fi
        if [   "$glo_module" == "Comm" ];      then glo_module="Comm_Mod";      fi
        if [   "$glo_module" == "Sort" ];      then glo_module="Sort_Mod";      fi
        if [   "$glo_module" == "Flow" ];      then glo_module="Field_Mod";     fi
        if [   "$glo_module" == "File" ];      then glo_module="File_Mod";      fi
        if [   "$glo_module" == "Profiler" ];  then glo_module="Profiler_Mod";  fi
        if [   "$glo_module" == "Convert" ];   then glo_module="Convert_Mod";   fi
        if [   "$glo_module" == "String" ];    then glo_module="String_Mod";    fi
      fi
      called_modules[$proc]=$glo_module
      called_procedures[$proc]=$glo_procedure
    done

    #------------------------------------------------------------------
    #   Print out the name of the module you are currently analysing
    #------------------------------------------------------------------
    if [ ! -z "$called_procedures" ]; then

      # It is in a module, print her LIGHT_GREEN
      if [ $module_in_which_you_seek ]; then
        echo "${indent}""-----------------------------------------------------------------------------------"
        echo -e ${LIGHT_GREEN}"${indent}"+ $procedure_name_you_seek "("$this_level")"${RESET}" calls: ""\t $module_in_which_you_seek"

      # If it is a global or external function, print it in light cyan
      else
        echo -e ${LIGHT_CYAN}"${indent}"+ $procedure_name_you_seek "("$this_level")"${RESET}" calls: ""\t $module_in_which_you_seek"
      fi
    else
      if [ $module_in_which_you_seek ]; then
        echo -e ${LIGHT_BLACK}"${indent}"⨯ $procedure_name_you_seek "("$this_level")"${RESET}"\t $module_in_which_you_seek"
      else
        echo -e ${LIGHT_RED}"${indent}"⨯ $procedure_name_you_seek "("$this_level")"${RESET}"\t $module_in_which_you_seek"
      fi
    fi

    # Increase indend for the next level by appending 6 spaces to it
    local indent="${indent}"'      '

    #-------------------------------------------------------------
    #   If the list of called procedures in not empty, carry on
    #-------------------------------------------------------------
    if [ ! -z "$called_procedures" ]; then

      # Print the procedures which are called from here
      for proc in "${!called_procedures[@]}"; do

        # This seems to be the only place from which you can print
        # non-member / external functions which don't make any calls.  Period!
        if [ ! ${called_modules[proc]} ]; then
          echo -e "${indent}"• ${LIGHT_CYAN}${called_procedures[proc]}${RESET}" \t (global | external | end of line)"  "mods: " ${called_modules[proc]}
        fi

        #-------------------------------------------------------#
        #                                                       #
        #   The very important recursive call to its own self   #
        #                                                       #
        #-------------------------------------------------------#
        extract_call_graph ${called_procedures[proc]} ${called_modules[proc]}
      done
    fi
  fi
}

#------------------------------------------------------------------------------
#  Three command line arguments are sent - process the second and pass them on
#------------------------------------------------------------------------------
if [ $3 ]; then
  if [ $2 == "-e" ]; then
    exclude_dir=""
    if [ "$3" ]; then
      exclude_dir="$3"
    fi
    extract_call_graph $1
  else
    print_usage
  fi

#---------------------------------------------------------------------
#  One command line argument is sent - must be the procedure name
#  Use the names without extension - say Grid_Mod, Convert_Mod ...
#---------------------------------------------------------------------
elif [ $1 ]; then
  extract_call_graph $1

#-----------------------------------------------------------------------
#  Wrong number of command line argument is sent - describe the usage
#-----------------------------------------------------------------------
else
  print_usage
fi

