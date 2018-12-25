import os
import sys
import getpass
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.utils import COMMASPACE, formatdate


def send_mail(send_from, passwd, send_to, mail_cc, subject, text, attachment=None, 
        smtp_server = 'smtp.126.com'):
    msg = MIMEMultipart()
    msg['Subject'] = subject
    msg['From'] = send_from
    msg['To'] = COMMASPACE.join(send_to)
    msg['Cc'] = COMMASPACE.join(mail_cc)

    msg.attach(MIMEText(text))

    # 添加附件
    if attachment is not None:
        with open(attachment, 'rb') as fh:
            part = MIMEApplication(
                    fh.read(),
                    Name=os.path.basename(attachment)
                    )
        part['Content-Disposition'] = 'attachment; filename=%s' % os.path.basename(attachment)
        msg.attach(part)

    try:
        server = smtplib.SMTP()
        server.set_debuglevel(1)
        server.connect(smtp_server)
        server.ehlo()
        server.starttls()
        server.ehlo()
        # login
        server.login(from_addr, passwd)
        # send mail
        server.sendmail(send_from, send_to + mail_cc, msg.as_string())
        server.close()
        print('邮件发送成功!')
    except Exception as e:
        print('邮件发送失败，', e)

if __name__ == '__main__':
    f = sys.argv[1]
    from_addr = 'wujunjames@126.com'
    passwd = getpass.getpass('Your email password:')
    to_addr = input('send to email:')
    mail_cc = input('Copy to email:')
    subject='引物合成'
    text = '你好，附件是需要合成的引物列表，请知悉。'
    send_mail(from_addr, passwd, [to_addr], [mail_cc], subject, text, attachment=f)


