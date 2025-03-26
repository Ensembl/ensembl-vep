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
# apt-get update
apt-get install -y mysql-server-5.7 mysql-client-5.7 mysql-client-core-5.7
mysql -e 'SET GLOBAL local_infile=1;'

rm mysql-apt-config_0.8.12-1_all.deb