## uninstall pre-installed mysql 8.0

sudo sudo apt-get remove --purge mysql*
sudo apt-get autoremove
sudo apt-get autoclean
sudo rm -rf /var/lib/mysql
sudo rm -rf /etc/mysql

## install mysql 5.7

# set deb-conf database answer so we are not prompted when configuring mysql-apt-conf
export DEBIAN_FRONTEND=noninteractive
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/repo-codename select bionic'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/repo-distro select ubuntu'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/repo-url string http://repo.mysql.com/apt/'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/select-preview select Disabled'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/select-product select Ok'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/select-server select mysql-5.7'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/select-tools select Disabled'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/unsupported-platform select abort'

# set root user password
sudo debconf-set-selections <<< "mysql-community-server mysql-community-server/root-pass password "
sudo debconf-set-selections <<< "mysql-community-server mysql-community-server/re-root-pass password "
sudo debconf-set-selections <<< "mysql-server-5.7 mysql-server/root_password password "
sudo debconf-set-selections <<< "mysql-server-5.7 mysql-server/root_password_again password "

# download and install mysql-apt-conf
wget http://dev.mysql.com/get/mysql-apt-config_0.8.12-1_all.deb
sudo dpkg -i mysql-apt-config_0.8.12-1_all.deb

# (NOT NEEDED CURRENTLY) install some dependencies that are not avaialable on focal repo
# wget http://archive.ubuntu.com/ubuntu/pool/main/liba/libaio/libaio1_0.3.110-5_amd64.deb
# sudo dpkg -i libaio1_0.3.110-5_amd64.deb
# wget http://archive.ubuntu.com/ubuntu/pool/universe/n/ncurses/libtinfo5_6.3-2_amd64.deb
# sudo dpkg -i libtinfo5_6.3-2_amd64.deb

# install public key otherwise will complain about mysql repo
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys B7B3B788A8D3785C
sudo apt-get update

# install mysql 5.7
sudo apt-cache policy mysql-server
sudo apt-get install -y mysql-server=5.7* mysql-client=5.7* libmysqlclient-dev

# start and configure mysql server
sudo apt-get install -y debconf-utils
sudo debconf-get-selections | grep mysql
sudo systemctl start mysql
sudo mysql -u root <<END
    CREATE USER 'travis'@'127.0.0.1' IDENTIFIED BY '';
    GRANT ALL PRIVILEGES ON *.* TO 'travis'@'127.0.0.1';
    FLUSH PRIVILEGES;
END

## clean up

rm mysql-apt-config_0.8.12-1_all.deb
# rm libaio1_0.3.110-5_amd64.deb
# rm libtinfo5_6.3-2_amd64.deb