"""Workspace locations and user-facing file extensions."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
PROJECTS = ROOT / 'projects'
PROCESSING = ROOT / 'processing'
WORK = ROOT / 'work'
EXAMPLES = ROOT / 'examples'
PROJECT_FILTER = 'Projekt DOSY Lab (*.proj)'
PROFILE_FILTER = 'Profil processingu DOSY Lab (*.list)'


def dialog_directory(directory):
    directory.mkdir(parents=True, exist_ok=True)
    return str(directory)


def require_extension(path, suffix):
    path = Path(path)
    if path.suffix.lower() != suffix:
        raise ValueError(f'Wybierz plik z rozszerzeniem {suffix}.')
    return path


def choose_save_file(parent, caption, directory, filename, file_filter):
    from PySide6.QtWidgets import QFileDialog

    dialog = QFileDialog(parent, caption, dialog_directory(directory), file_filter)
    dialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
    dialog.setDefaultSuffix(Path(filename).suffix[1:])
    dialog.selectFile(filename)
    return dialog.selectedFiles()[0] if dialog.exec() else None
