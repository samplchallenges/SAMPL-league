# Generated by Django 3.2 on 2021-04-20 21:09

from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('core', '0004_blobvalue_floatvalue_prediction_textvalue'),
    ]

    operations = [
        migrations.CreateModel(
            name='Evaluation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('exit_status', models.IntegerField()),
                ('log_stdout', models.TextField(blank=True)),
                ('log_stderr', models.TextField(blank=True)),
                ('input_element', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='core.inputelement')),
            ],
        ),
        migrations.CreateModel(
            name='SubmissionRun',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('digest', models.CharField(max_length=255)),
                ('data_privacy_level', models.CharField(choices=[('PUBLIC', 'Public'), ('PRIVATE', 'Private')], default='PRIVATE', max_length=25)),
                ('status', models.CharField(choices=[('FAILURE', 'Failure'), ('SUCCESS', 'Success')], max_length=25)),
                ('submission', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='core.submission')),
            ],
            options={
                'get_latest_by': ('created_at',),
                'abstract': False,
            },
        ),
        migrations.RemoveField(
            model_name='submissionevaluation',
            name='submission',
        ),
        migrations.DeleteModel(
            name='SubmissionResult',
        ),
        migrations.AddField(
            model_name='evaluation',
            name='submission_run',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='core.submissionrun'),
        ),
        migrations.AddField(
            model_name='prediction',
            name='evaluation',
            field=models.ForeignKey(default=1, on_delete=django.db.models.deletion.CASCADE, to='core.evaluation'),
            preserve_default=False,
        ),
        migrations.AlterUniqueTogether(
            name='prediction',
            unique_together={('evaluation', 'key')},
        ),
        migrations.AlterUniqueTogether(
            name='evaluation',
            unique_together={('submission_run', 'input_element')},
        ),
        migrations.RemoveField(
            model_name='prediction',
            name='input_element',
        ),
        migrations.RemoveField(
            model_name='prediction',
            name='submission_evaluation',
        ),
        migrations.DeleteModel(
            name='SubmissionEvaluation',
        ),
    ]