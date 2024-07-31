## 1.简介
Git 是一个开源的分布式版本控制系统，用于敏捷高效地处理任何或小或大的项目。使用Git进行版本控制，可以方便的查看修改存储。

### 2.下载git
* 打开[Git官网](https://git-scm.com/downloads)下载。
* 安装Git.
### 3.远程仓库配置
远程仓库有很多，比如github，gitlab，国内的gitee，局域网自建git服务器等。这里推荐使用GitHub.com，作为代码远程托管仓库。

* 使用github.com网站，注册GitHub账号（若打不开网址，建议开启科学上网使用代理。）
* 在自己的账号下新建一个仓库，名称任意这里为了演示命名为Test_Repo。
### 4.使用SSH生成配置连接本地到远端
如果使用https方式在每次git pull和git push时需要输入用户名和密码来验证登录，操作繁琐且麻烦。使用SSH方式可无需验证登录，操作简便，推荐使用SSH的方式使本地仓库与远程仓库建立链接。

* 在本地任意目录下，右键鼠标选择<span style="color:green;">Open Git Bash here</span>，可以看到打开了一个Bash终端。输入以下命令：
```bash
cd ~ #进入到个人用户目录
git config --global user.name "注册名" #任意名称，使用英文
git config --global user.email "注册邮箱" #任意邮箱，使用英文
ssh-keygen -t rsa -C "自己的邮箱" #生成ssh公钥，后续需要按三次回车键
```
连按三次 Enter，最后会生成两个文件：id_rsa 和 id_rsa.pub。SSH文件存放在C:/User/用户/.ssh下，id_rsa为私钥，不能泄露，id_rsa.pub为公钥，用于建立链接。
* github.com填写SSH Key
> 具体为点击个人账户头像，找到settings，接着点击左侧SSH and GPG keys，点击New SSH key。title可以任意，key填写id_rsa.pub的全部内容。点击Add SSH key。
* 验证是否连接成功
在bash终端中输入
```bash
ssh -T git@github.com
```
注意：首次设置需要输入yes，查看是有successfully authenticated字样，有就表明关联成功。
### 5.将本地仓库推送到远端仓库
在如下目录D:\DEVELOP\下新建文件夹命名为firstTry。在firstTry中右键Open Git Bash here，打开一个终端输入：
```bash
git init #初始化本地仓库。
```
这个命令初始化当前文件夹为本地的一个仓库，这时会在当前目录下产生一个.git的文件夹，就是用来跟踪和管理版本库的。接下来去自己的GitHub账号下，将刚才新建的远程仓库Test_Repo的链接复制下来：
```bash
git remote add origin git@github.com:rongguangyiii/Test_Repo.git
```
这样就建立了本地仓库和远端仓库的关联。
### 6.使用git命令
建立了本地仓库和远程仓库的连接，就可以使用git的命令将本地更改推送到远端。首次推送仓库里的内容使用：
```bash
git status #查看当前分支，本地仓库的状态
git add . #将更改添加到暂存库
git commit -m "First commit1" #将更改提交到本地仓库，并添加说明
git status 
git push -u origin master #将本次修改提交到远端
git log #查看提交log及分支提交情况
```
首次需要加<span style="color:red;">-u</span>这个参数。当远端有了仓库以后，就不需要加-u了。
### 7. git常用命令学习
* 详细学习教程参考Git官网的文档(https://git-scm.com/book/zh/v2)
* [git命令速查表](https://cheatsheet.wang/)

Author : Guangying

Date :  2024.04.07
