## 1.���
Git ��һ����Դ�ķֲ�ʽ�汾����ϵͳ���������ݸ�Ч�ش����κλ�С������Ŀ��ʹ��Git���а汾���ƣ����Է���Ĳ鿴�޸Ĵ洢��

### 2.����git
* ��[Git����](https://git-scm.com/downloads)���ء�
* ��װGit.
### 3.Զ�ֿ̲�����
Զ�ֿ̲��кܶ࣬����github��gitlab�����ڵ�gitee���������Խ�git�������ȡ������Ƽ�ʹ��GitHub.com����Ϊ����Զ���йֿܲ⡣

* ʹ��github.com��վ��ע��GitHub�˺ţ����򲻿���ַ�����鿪����ѧ����ʹ�ô�����
* ���Լ����˺����½�һ���ֿ⣬������������Ϊ����ʾ����ΪTest_Repo��
### 4.ʹ��SSH�����������ӱ��ص�Զ��
���ʹ��https��ʽ��ÿ��git pull��git pushʱ��Ҫ�����û�������������֤��¼�������������鷳��ʹ��SSH��ʽ��������֤��¼��������㣬�Ƽ�ʹ��SSH�ķ�ʽʹ���زֿ���Զ�ֿ̲⽨�����ӡ�

* �ڱ�������Ŀ¼�£��Ҽ����ѡ��<span style="color:green;">Open Git Bash here</span>�����Կ�������һ��Bash�նˡ������������
```bash
cd ~ #���뵽�����û�Ŀ¼
git config --global user.name "ע����" #�������ƣ�ʹ��Ӣ��
git config --global user.email "ע������" #�������䣬ʹ��Ӣ��
ssh-keygen -t rsa -C "�Լ�������" #����ssh��Կ��������Ҫ�����λس���
```
�������� Enter���������������ļ���id_rsa �� id_rsa.pub��SSH�ļ������C:/User/�û�/.ssh�£�id_rsaΪ˽Կ������й¶��id_rsa.pubΪ��Կ�����ڽ������ӡ�
* github.com��дSSH Key
> ����Ϊ��������˻�ͷ���ҵ�settings�����ŵ�����SSH and GPG keys�����New SSH key��title�������⣬key��дid_rsa.pub��ȫ�����ݡ����Add SSH key��
* ��֤�Ƿ����ӳɹ�
��bash�ն�������
```bash
ssh -T git@github.com
```
ע�⣺�״�������Ҫ����yes���鿴����successfully authenticated�������оͱ��������ɹ���
### 5.�����زֿ����͵�Զ�˲ֿ�
������Ŀ¼D:\DEVELOP\���½��ļ�������ΪfirstTry����firstTry���Ҽ�Open Git Bash here����һ���ն����룺
```bash
git init #��ʼ�����زֿ⡣
```
��������ʼ����ǰ�ļ���Ϊ���ص�һ���ֿ⣬��ʱ���ڵ�ǰĿ¼�²���һ��.git���ļ��У������������ٺ͹���汾��ġ�������ȥ�Լ���GitHub�˺��£����ղ��½���Զ�ֿ̲�Test_Repo�����Ӹ���������
```bash
git remote add origin git@github.com:rongguangyiii/Test_Repo.git
```
�����ͽ����˱��زֿ��Զ�˲ֿ�Ĺ�����
### 6.ʹ��git����
�����˱��زֿ��Զ�ֿ̲�����ӣ��Ϳ���ʹ��git��������ظ������͵�Զ�ˡ��״����Ͳֿ��������ʹ�ã�
```bash
git status #�鿴��ǰ��֧�����زֿ��״̬
git add . #��������ӵ��ݴ��
git commit -m "First commit1" #�������ύ�����زֿ⣬�����˵��
git status 
git push -u origin master #�������޸��ύ��Զ��
git log #�鿴�ύlog����֧�ύ���
```
�״���Ҫ��<span style="color:red;">-u</span>�����������Զ�����˲ֿ��Ժ󣬾Ͳ���Ҫ��-u�ˡ�
### 7. git��������ѧϰ
* ��ϸѧϰ�̳̲ο�Git�������ĵ�(https://git-scm.com/book/zh/v2)
* [git�����ٲ��](https://cheatsheet.wang/)

Author : Guangying

Date :  2024.04.07
