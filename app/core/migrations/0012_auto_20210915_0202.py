# Generated by Django 3.2.7 on 2021-09-15 02:02

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0011_submissionarg'),
    ]

    operations = [
        migrations.AlterField(
            model_name='evaluation',
            name='status',
            field=models.CharField(choices=[('FAILURE', 'Failure'), ('SUCCESS', 'Success'), ('PENDING', 'Pending'), ('RUNNING', 'Running'), ('CANCELLED', 'Cancelled')], default='PENDING', max_length=25),
        ),
        migrations.AlterField(
            model_name='submissionrun',
            name='status',
            field=models.CharField(choices=[('FAILURE', 'Failure'), ('SUCCESS', 'Success'), ('PENDING', 'Pending'), ('RUNNING', 'Running'), ('CANCELLED', 'Cancelled')], max_length=25),
        ),
    ]
