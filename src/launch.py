"""
One-stop script to set up and launch DVR process on multiple VMs!

Pre-requisites:
- Install paramiko and scp:
   pip install paramiko
   pip install scp
- Make sure SSH keys are added to each VM.
- SSH into each VM at least once prior to running this script so that the IPs will be added to the known_hosts file.
  If you don't do this, you will get an error saying <IP> not found in known_hosts

Notes:
- If you get the error "Host key for server '34.136.55.60' does not match: got ..., expected ...
  then you need to remove the entry for that IP in the known_hosts file (simply delete the line(s) corresponding to that IP)
  and re-add it to known_hosts by SSHing into the VM with that IP.
"""
from paramiko import SSHClient
from scp import SCPClient
import threading
from queue import Queue, Empty
import sys
import time

# global variables
uni = "rdn2108"
key_filename = "/Users/rdn/.ssh/csee4119-s24"
vm_ex_ips = [
'34.134.176.42',
'35.192.179.63',
'35.225.219.240'
]
internal_port = 50000
external_port = 60000
topology_file = 'topology.dat'


class bcolors
    """
    A utility class to define colors for printing to the console

    Source: https://stackoverflow.com/questions/287871/how-do-i-print-colored-text-to-the-terminal
    """
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    RESET = '\x1b[0m'


def kill_overlays():
   """
   Kill any overlay processes running on the VMs
   """
   try:
      ssh = SSHClient()
      ssh.load_system_host_keys()

      print(bcolors.OKBLUE + "Killing any running overlay processes on VMs" + bcolors.RESET)
      for vm_ex_ip in vm_ex_ips:
         ssh.connect(vm_ex_ip, username=uni, key_filename=key_filename)
         kill_stdin, kill_stdout, kill_stderr = ssh.exec_command("pkill -9 overlay")
         ssh.close()
   except KeyboardInterrupt:
      print(bcolors.OKBLUE + "Bye!" + bcolors.RESET)
      sys.exit()

def load_files():
   """
   SCPs the overlay exectuable, dvr.py and topology dat to all VMs
   """
   try:
      ssh = SSHClient()
      ssh.load_system_host_keys()

      print(bcolors.OKBLUE + "SCPing overlay, dvr.py, and topology.dat to VMs" + bcolors.RESET)
      for vm_ex_ip in vm_ex_ips:
         ssh.connect(vm_ex_ip, username=uni, key_filename=key_filename)
         with SCPClient(ssh.get_transport()) as scp:
            scp.put('overlay', 'overlay')
            scp.put('dvr.py', 'dvr.py')
            scp.put('topology.dat', 'topology.dat')
            stdin, stdout, stderr = ssh.exec_command("chmod +x overlay")
         ssh.close()
   except KeyboardInterrupt:
      print(bcolors.OKBLUE + "Bye!" + bcolors.RESET)
      sys.exit()

def start_overlay(vm_name, vm_ex_ip, ready_queue):
   """
   Start overlay process in a given VM. Prints any output of overlay process running on the VM to console.

   Args:
      vm_name : A string or int representing the name of the VM, only used for logging.
      vm_ex_ip : A string representing the external IP of the VM
      ready_queue : A queue object for start_overlay to signal to the main thread
                  when the overlay process is ready to receive network connections
   """

   ssh = SSHClient()
   ssh.load_system_host_keys()
   ssh.connect(vm_ex_ip, username=uni, key_filename=key_filename)
   o_run_stdin, o_run_stdout, o_run_stderr = ssh.exec_command( f"stdbuf -oL ./overlay {external_port} {internal_port} {topology_file} 2> \&1 &") # start overlay process
   time.sleep(2) # wait a bit for overlay to start, then ensure that is started successfully, otherwise exit because won't ever work
   _, overlay_process_stoud, _ = ssh.exec_command("pidof overlay")
   if len(overlay_process_stoud.read()) == 0:
      print(bcolors.FAIL + f"VM  {vm_name} [Overlay]: Couldn't start the overlay on this VM. Please temrminate and rerun the program." + bcolors.RESET)
   for line in o_run_stdout:
      if "Overlay: waiting for connection from the network process..." in line:
         ready_queue.put("READY")
      print(f"VM  {vm_name} [Overlay]: {line.rstrip()}")
   ssh.close()


def start_dvr(vm_name, vm_ex_ip):
   """
   Start the DVR process in a given VM. Prints any output of dvr.py process running on the VM to console.

   Args:
      vm_name : A string or int representing the name of the VM, only used for logging.
      vm_ex_ip : A string representing the external IP of the VM
   """
   ssh = SSHClient()
   ssh.load_system_host_keys()
   ssh.connect(vm_ex_ip, username=uni, key_filename=key_filename)
   dvr_run_stdin, dvr_run_stdout, dvr_run_stderr = ssh.exec_command( f"stdbuf -oL python3 dvr.py {external_port}")
   for line in dvr_run_stdout:
      print(f"VM  {vm_name} [DVR] {line.rstrip()}")
   ssh.close()


def launch():
   """
   Start overlay and DVR processes in all VMs. Close gracefully on KeyboardInterrupt.
   """
   try:
      threads = []

      # start overlay processes in each machine
      overlay_ready = [0 for i in range(len(vm_ex_ips))]
      ready_queues = [Queue() for i in range(len(vm_ex_ips))]
      for vm_name, vm_ex_ip in enumerate(vm_ex_ips):
         t = threading.Thread(target = start_overlay, args=(vm_name, vm_ex_ip, ready_queues[vm_name], ))
         t.daemon = True
         t.start()
         threads.append(t)

      # wait until each overlay process is ready to receive network connections
      while not all(overlay_ready):

         for i, q in enumerate(ready_queues):
            try:
               if q.get(True, 0.1) == "READY":
                  overlay_ready[i] = 1
                  print(bcolors.OKGREEN + f"VM {i} ready for connections." + bcolors.RESET)
            except Empty:
               pass


      # start DVR processes in each machine
      print(bcolors.OKBLUE + "Starting DVR scripts" + bcolors.RESET)
      for vm_name, vm_ex_ip in enumerate(vm_ex_ips):
         t = threading.Thread(target = start_dvr, args=(vm_name, vm_ex_ip, ))
         t.start()
         threads.append(t)

   except KeyboardInterrupt:
      print(bcolors.OKBLUE + "Bye!" + bcolors.RESET)
      sys.exit()


if __name__ == "__main__":
   kill_overlays() # you will always need to run this before load_files() or launch() as a precaution
   #load_files() # comment this out if you have already SCP'd the files
   launch()
