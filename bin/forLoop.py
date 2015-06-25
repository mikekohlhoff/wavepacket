from numpy import arange
from string import Template, split

import os
import subprocess
import sys
import shutil
import threading
import time
import math
import time

# Loop through a range of one variable (likely to be velocity), creating an input file for each and
# execute CWDVR. The script limits the number of concurrently executing CWDVR
# simulations to the value in max_threads.

# conversion vel in atomic units to m/s
v_hatree = 2.1876912633E6
v_step = 400
v_vals = (arange(v_step,3200+1,v_step)*-1*math.sin(20)/v_hatree)
print(arange(v_step,3200+1,v_step))
time.sleep(3)
v_vals = (('%e' % i).replace('e', 'd') for i in v_vals)

# Base path is the directory in which each sub-directory is created. Run path
# is a template for the name of each working directory. Each value of i and j
# are used to create the directory name through substitution of ${i} and ${j}
# placeholders respectively. The template file is the source for the input file
# that will be copied to each working directory. The markers ${i}
# replaced with each value from the loop.
base_path = r'../dataout/n=5k=-4'
run_path  = r'n=5k=-4v=${v}'
config_file = r'config.info'

# CWDVR execution. Specify the path to the cwdvr executable, and specify the
# maximum number of copies to run simultaneously.
cwdvr_command = os.path.abspath(r'./cwdvr')
# ccmd_command = os.path.expanduser(ccmd_command)
# TPS7/8 have 8 cores with max 2 threads
max_threads = 4

#======================================================================

# Read base input file and convert to a Python template.
try:
    basefile = open(config_file)
    params_template = Template(basefile.read())
    path_template = Template(run_path)
    basefile.close()
except IOError as e:
    print "Error opening template file %s\nerror(%s) : %s" \
            % (template_file, e[0], e[1])
    sys.exit()


def do_loop():
    """
    Perform the loop over parameters. During each iteration, the name of the
    working directory is built from the current values of v and the
    template stored in run_path. This directory is created if it doesn't
    already exist. Placeholders marked ${v} CWDVR input file
    template_file are replaced with the current values and the file is written
    to the new working drectory, overwriting one already there.

    Each working directory is passed to a ThreadCounter object that handles
    calling the cwdvr executable and limiting the number of concurrently running
    copies.
    """

    for v_val in v_vals:
        # Build a dictionary of replacements for the placeholders
        replacements = dict(v=v_val)
        # Start a thread to execute CWDVR when ready.
        exe = ThreadCounter(replacements)
        exe.start()

#===============================================================================

# Create a semaphore to keep count of the number of running copies.
thread_limit = threading.BoundedSemaphore(max_threads)

class ThreadCounter(threading.Thread):
    """
    Class to handle executing a new copy of CWDVR when ready. The execution is
    blocked until the semaphore can be acquired, then the command-line for
    CWDVR is passed to a subprocess for execution.
    """
    def __init__(self, replacements):
        threading.Thread.__init__(self)
        self.replacements = replacements
        print self.replacements
        # convert to fortran double precision real
        # self.replacements = repr(self.replacements)

    def run(self):
        thread_limit.acquire()
        try:
            # Substitute placeholders in file path and input file templates.
            # Fail and quit the script if we are left with anything
            # un-substituted.
            new_params = params_template.substitute(self.replacements)
            new_path   = path_template.substitute(self.replacements)


            # Concatentate the base path with the new path name and create this
            # directory if it doesn't exist. Fail and quit the script if this
            # causes problems.
            target_path = os.path.join(base_path, new_path)
            print target_path
            if not os.path.exists(target_path):
                os.makedirs(target_path)

            # Write the new, substituted, input file to the working directory
            newfile = open(os.path.join(target_path, "config.info"), 'w')
            newfile.write(new_params)
            newfile.close()


            print "\n********************"
            print "Thread starting in %s" % new_path
            print "********************\n"
            cmd = split(cwdvr_command)
            proc = subprocess.Popen(cmd, cwd=target_path)
            proc.communicate()

        except KeyError as ke:
            print "Unexpected key in file"
            print str(ke)
            sys.exit(1)
        except os.error as e:
            print "Can't create new working directory. Exiting."
            print e.strerror
            sys.exit(1)
        except IOError as e:
            print "Error writing input file %s\nerror(%s) : %s" \
                    % (target_path, e.errno, e.errstr)
        finally:
            thread_limit.release()

def main():
    """
    Call the primary function when the script is run from the command line.
    """
    do_loop()

if __name__ == "__main__":
    main()
