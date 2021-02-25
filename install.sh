#!/bin/bash

# check if linux/Mac or other
uname_out="$(uname -s)"
case $uname_out in
	Linux*)	machine=Linux;;
	Darwin*)	machine=Mac;;
	*)		machine=$uname_out;echo "unknown machine type: $uname_out, this installer probably won't work...";
esac

cp progs/ANSURR .

# install for all users?
#if [ "$(whoami)" != "root" ]; then
while true; do
	read -p " -> would you like to install ANSURR for all users? [recommended - requires root privileges] (y/n) " yn
	case $yn in
	[Yy]* ) sudo true; local_install=0; break;;
	[Nn]* ) local_install=1; break;;
	* ) echo "please answer yes or no";;
	esac
done
#fi

# set which version of python ANSURR should use
python_version=python
echo -n 'default version of python on this machine is '; python -V

while true; do
    read -p " -> would you like ANSURR to use this version? [python3.x is recommended] (y/n) " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) echo -n 'other versions of python installed on this machine include: '
		echo -n $(ls /usr/bin | grep ^python | grep -v config)' '; echo $(ls /usr/local/bin | grep ^python | grep -v config)
		echo -n ' -> please type which to use from the list above (or type the path if not in the list): '
		while true; do
			read python_version
			if [ -z $(command -v $python_version) ]; then
				echo "$python_version not recognised" 
				echo -n " -> please try again: "
			else
				break
			fi
		done
		break;;
        * ) echo "please answer yes or no";;
    esac
done

# check python version for required modules 
echo "checking $python_version has required modules (numpy, scipy, matplotlib)"; 
missing_modules=$($python_version progs/check_python.py)
if [ ! -z "$missing_modules" ]; then
	if [ $local_install == "1" ]; then
		echo "could not find$missing_modules, you could try using pip to install by running: $python_version -m pip install$missing_modules --user"
	else
		echo "could not find$missing_modules, you could try using pip to install by running: $python_version -m pip install$missing_modules"
	fi
	while true; do
	    read -p " -> would you like to try this now? (y/n) " yn
	    case $yn in
		[Yy]* ) if [ -z $(command -v pip) ]; then
				echo "ERROR pip doesn't seem to be installed, install$missing_modules and re-run install.sh"; exit 1
			else
				if [ $local_install == "1" ]; then
					$python_version -m pip install$missing_modules --user
				else
					$python_version -m pip install$missing_modules
				fi
				missing_modules=$($python_version progs/check_python.py)
				if [ -z "$missing_modules" ]; then
					echo " -> found required modules"
				else
					echo "ERROR could not find$missing_modules, install$missing_modules and re-run install.sh"; exit 1 
				fi
					
			fi
			break;;
		[Nn]* ) echo "ERROR install$missing_modules and re-run install.sh"; exit 1; break;;
		* ) echo "please answer yes or no";;
	    esac
	done
else
	echo " -> found required modules"
fi

# check java is installed
echo "ANSURR can optionally use PANAV to re-reference chemical shifts, which requires Java, checking to see if Java is installed"
if [ -z $(command -v java) ]; then
	while true; do
		if [ $machine == "Mac" ]; then
	    	echo " -> Java not found, you could try to install Java by running: sudo brew install default-jre"
			java_cmd='sudo brew install default-jre'
		else
			echo " -> Java not found, you could try to install Java by running: sudo apt install default-jre"
			java_cmd='sudo apt install default-jre'
		fi
	    read -p " -> would you like to try this now? (y/n) " yn
	    case $yn in
		[Yy]* ) $java_cmd; 
			if [ -z $(command -v java) ]; then
				echo " -> Java failed to install"
			else
				echo " -> Java installed"
			fi
			break;;
		[Nn]* ) break;;
		* ) echo "please answer yes or no";;
	    esac
	done
else
	echo " -> found Java" 	
fi

# install
#sed -i '/local_install=/d' ANSURR
if [ $machine == "Mac" ]; then
	sed -i '' '/ANSURR_DIR=/d' ANSURR
	sed -i '' '/python_version=/d' ANSURR
else
	sed -i '/ANSURR_DIR=/d' ANSURR
	sed -i '/python_version=/d' ANSURR
fi

if [ $local_install == "1" ]; then
    ANSURR_DIR=$(pwd)
    if [ $machine == "Mac" ]; then
		sed "1a\\
                python_version=$python_version\\
                " ANSURR > ANSURR_tmp1
                sed "2a\\
                ANSURR_DIR=$ANSURR_DIR\\
                " ANSURR_tmp1 > ANSURR_tmp
                rm ANSURR_tmp1
	else
		sed "1 a\python_version=$python_version\nANSURR_DIR=$ANSURR_DIR" ANSURR > ANSURR_tmp
	fi
	chmod +x ANSURR_tmp
	mv ANSURR_tmp ansurr_local
	rm ANSURR
	if [ -z $(command -v ./ansurr_local) ]; then
		echo "ANSURR failed to install" 
		echo
	else
		echo "ANSURR installed successfully"
		echo 
		./ansurr_local -h
		echo
	fi
else
    ANSURR_DIR='/usr/local/lib/ansurr'
    if [ $machine == "Mac" ]; then
		sed "1a\\
		python_version=$python_version\\
		" ANSURR > ANSURR_tmp1
		sed "2a\\
                ANSURR_DIR=$ANSURR_DIR\\
                " ANSURR_tmp1 > ANSURR_tmp
		rm ANSURR_tmp1
	else
		sed "1 a\python_version=$python_version\nANSURR_DIR=$ANSURR_DIR" ANSURR > ANSURR_tmp
	fi
	chmod +x ANSURR_tmp
	sudo mkdir -p /usr/local/lib/ansurr
	if [ ! -d /usr/local/bin ]; then
        	sudo mkdir -p /usr/local/bin
	fi
	sudo cp -r progs /usr/local/lib/ansurr/
	sudo cp -r lib /usr/local/lib/ansurr/
	sudo mv ANSURR_tmp /usr/local/bin/ansurr
	rm ANSURR
	if [ -z $(command -v ansurr) ]; then
		echo "ANSURR failed to install" 
		echo
	else
		echo "ANSURR installed successfully" 
		echo
		ansurr -h
		echo
	fi

fi



