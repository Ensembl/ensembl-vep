# uninstall pre-installed mysql 8.0
sudo sudo apt-get remove --purge mysql*
sudo apt-get autoremove
sudo apt-get autoclean
sudo rm -rf /var/lib/mysql
sudo rm -rf /etc/mysql

# install mysql 5.7
export DEBIAN_FRONTEND=noninteractive
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/repo-codename select focal'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/repo-distro select ubuntu'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/repo-url string http://repo.mysql.com/apt/'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/select-preview select '
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/select-product select Ok'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/select-server select mysql-5.7'
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/select-tools select '
sudo debconf-set-selections <<< 'mysql-apt-config mysql-apt-config/unsupported-platform select abort'

wget http://dev.mysql.com/get/mysql-apt-config_0.8.12-1_all.deb
sudo dpkg -i mysql-apt-config_0.8.12-1_all.deb
sudo apt-get update
sudo apt-get install -y mysql-server=5.7* mysql-client=5.7*
sudo systemctl start mysql
sudo mysql -u root -e "CREATE USER 'travis'@'127.0.0.1' IDENTIFIED BY ''; FLUSH PRIVILEGES;"
mysql -e 'SET GLOBAL local_infile=1;'

rm mysql-apt-config_0.8.12-1_all.deb