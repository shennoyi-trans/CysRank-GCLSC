import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tools import launch


class LauncherTests(unittest.TestCase):
    def test_missing_and_wrong_versions(self):
        def version(name):
            if name == 'absent':
                raise launch.importlib.metadata.PackageNotFoundError(name)
            return {'correct': '1', 'wrong': '0'}[name]
        with patch.object(launch, 'required_versions', return_value={'correct': '1', 'wrong': '1', 'absent': '1'}), \
                patch.object(launch.importlib.metadata, 'version', side_effect=version):
            self.assertEqual(launch.missing_requirements(), ['wrong==1', 'absent==1'])

    def test_complete_install_does_not_reinstall(self):
        with patch.object(launch, 'missing_requirements', return_value=[]), patch.object(launch.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run:
            launch.install_dependencies()
        commands = [call.args[0] for call in run.call_args_list]
        self.assertFalse(any('install' in command for command in commands))
        self.assertTrue(any('check' in command for command in commands))
        self.assertTrue(any('--verify-imports' in command for command in commands))

    def test_install_uses_same_python_and_single_requirements(self):
        with patch.object(launch, 'missing_requirements', return_value=['biopython==1.88']), patch.object(launch.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run:
            launch.install_dependencies()
        commands = [call.args[0] for call in run.call_args_list]
        self.assertTrue(all(command[0] == launch.sys.executable for command in commands))
        self.assertEqual(next(command for command in commands if 'install' in command)[-2:], ['-r', str(launch.REQUIREMENTS)])

    def test_missing_transitive_dependency_triggers_repair(self):
        with patch.object(launch, 'missing_requirements', return_value=[]), patch.object(launch.subprocess, 'run', return_value=SimpleNamespace(returncode=1, stdout='missing transitive dependency')) as run:
            launch.install_dependencies()
        self.assertTrue(any('install' in call.args[0] for call in run.call_args_list))

    def test_external_environment_is_not_modified(self):
        with tempfile.TemporaryDirectory(prefix='project with spaces ') as directory, \
                patch.object(launch, 'ROOT', Path(directory)), patch.object(launch.sys, 'prefix', 'external'), \
                patch.object(launch.sys, 'base_prefix', 'external'), patch.object(launch.subprocess, 'run') as run:
            python = launch.managed_python()
            self.assertIn('--system-site-packages', run.call_args.args[0])
            self.assertEqual(python, Path(directory) / '.runtime/venv-313/Scripts/python.exe')

    def test_project_environment_is_reused(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(launch, 'ROOT', Path(directory)), \
                patch.object(launch.sys, 'prefix', str(Path(directory) / '.venv')), \
                patch.object(launch.sys, 'base_prefix', 'external'), patch.object(launch.subprocess, 'run') as run:
            self.assertEqual(launch.managed_python(), Path(launch.sys.executable))
            run.assert_not_called()

    def test_broken_environment_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(launch, 'ROOT', Path(directory)), \
                patch.object(launch.sys, 'prefix', 'external'), patch.object(launch.sys, 'base_prefix', 'external'), \
                patch.object(launch.subprocess, 'run'):
            previous = Path(directory) / '.runtime/venv-313'
            previous.mkdir(parents=True)
            (previous / 'keep.txt').write_text('preserve')
            self.assertNotEqual(launch.managed_python().parent.parent, previous)
            self.assertEqual((previous / 'keep.txt').read_text(), 'preserve')

    def test_missing_asset_fails_early(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(launch, 'ROOT', Path(directory)):
            with self.assertRaisesRegex(RuntimeError, 'app.py'):
                launch.verify_assets()

    def test_model_cache_reuse_download_and_offline_failure(self):
        from huggingface_hub.errors import LocalEntryNotFoundError
        for mode in ('project', 'global', 'download', 'offline'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                for name in ('config.json', 'pytorch_model.bin'):
                    (folder / name).write_text('fixture')
                def download(repo, name, **kwargs):
                    if mode == 'offline' or (mode == 'download' and kwargs.get('local_files_only', False)):
                        raise LocalEntryNotFoundError('fixture cache miss')
                    if mode == 'global' and kwargs['cache_dir'] is not None:
                        raise LocalEntryNotFoundError('fixture cache miss')
                    return str(folder / name)
                with patch('huggingface_hub.hf_hub_download', side_effect=download) as hub:
                    if mode == 'offline':
                        with self.assertRaisesRegex(RuntimeError, '缓存'):
                            launch.prepare_model(offline=True)
                    else:
                        self.assertEqual(launch.prepare_model(), folder)
                    online_calls = [call for call in hub.call_args_list if not call.kwargs.get('local_files_only', False)]
                    self.assertEqual(len(online_calls), 2 if mode == 'download' else 0)


if __name__ == '__main__':
    unittest.main()
