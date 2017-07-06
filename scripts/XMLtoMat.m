dir = "HRpQCT/";
output_dir = "HRpQCT/mats/";
filelist = readdir (dir)
for ii = 1:numel(filelist)
  ## skip special files . and ..
  if (regexp (filelist{ii}, "^\\.\\.?$") || regexp (filelist{ii}, ".m$"))
    continue;
  endif

  actual_file = filelist{ii};
  if(regexp (filelist{ii}, "Slices.xml$"))
    S = OpenXML(strcat(dir, actual_file));
    [~, pure_filename, ~] = fileparts(actual_file);
    filename = strcat(output_dir, pure_filename, ".mat");
    save("-mat-binary", filename, "S");
  else
    if(regexp (filelist{ii}, "Mask.xml$"))
      M = OpenMaskXML(strcat(dir, actual_file));
      [~, pure_filename, ~] = fileparts(actual_file);
      filename = strcat(output_dir, pure_filename, ".mat");
      save("-mat-binary", filename, "M");
    endif
  endif

endfor
