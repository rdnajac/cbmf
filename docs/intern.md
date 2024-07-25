# AWS guest access

| Username | Password | Console sign-in URL |
|-----------|----------|----------------------|
| guest     | ############### | [AWS Console](https://ferrando-lab.signin.aws.amazon.com/console) |

## Connect to EC2 instance

Go to the [EC2 Console](https://us-east-1.console.aws.amazon.com/ec2/home?region=us-east-1#Instances:instanceState=running)
and select the instance you want to connect to.  Click on the `Connect`,
then `EC2 serial console` tab, and `Connect` again.

![EC2 serial console](/images/ec2-serial-console.png)

If everything is working correctly, you should see a prompt like this:

![EC2 serial console login](/images/ec2-console-prompt.png)

> [!NOTE]
> If you encounter a login prompt, use the following credentials:
>
> - Username: <yourfirstname>
> - Password: <same as the AWS guest password>

