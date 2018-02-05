## Nimbus login
  ssh -i ~/Downloads/siobhon_egan_kp.pem ubuntu@146.118.65.16

Check for updates
  sudo apt-get update

*Adding a volume*
[Pawsey info](https://support.pawsey.org.au/documentation/display/US/Nimbus+-+Creating+and+Attaching+Volumes)

  root@test-instance:~# sudo fdisk -l /dev/vdc
  root@test-instance:~# sudo mkfs.ext4 /dev/vdc
  root@test-instance:~# sudo mkdir /data

Once volume has been added need to change permissions to allow you to use volume
  root@test-instance:~# sudo chown ubuntu:ubuntu /volumename

nohup ./script &

nohup ./script &> output.txt &

top #this will show what is currently running

kill #jobID #this wil kill the job
