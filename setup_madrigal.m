% download patched madrigal API, extracting under this directory.

cwd = fileparts(mfilename("fullpath"));
mad_dir = fullfile(cwd, "madmatlab");

if isfile(fullfile(mad_dir, "globalIsprint.m"))
  disp(mad_dir)
  addpath(mad_dir)
  return
end

error("please type in Terminal (not Matlab)\ngit -C %s submodule init --update", cwd)
