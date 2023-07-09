% download patched madrigal API, extracting under this directory.

cwd = fileparts(mfilename("fullpath"));
mad_dir = fullfile(cwd, "matlab-madrigal-main");

mad_archive = "remoteMatlabAPI.zip";
url = "https://github.com/semeter-lab/matlab-madrigal/archive/refs/heads/main.zip";

if isfile(fullfile(mad_dir, "globalIsprint.m"))
  disp(mad_dir)
  addpath(mad_dir)
  return
end

if ~isfile(mad_archive)
  disp(url + " => " + mad_archive)
  websave(mad_archive, url);
end

disp(mad_archive + " => " + mad_dir)
unzip(mad_archive, cwd)

addpath(mad_dir)
disp(mad_dir)
