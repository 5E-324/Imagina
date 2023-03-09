#pragma once

enum class FileType {
	Imagina,
	ImaginaText,
	Kfr,
	Kfp,
	Other,
};

void OpenFile(wchar_t *FileName, size_t ExtensionOffset);
void SaveFile(wchar_t *FileName, FileType Type);
void SaveImage(wchar_t *FileName);
void SaveRawPixelData(wchar_t *FileName);
