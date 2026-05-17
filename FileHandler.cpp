#include "FileHandler.H"

#include <algorithm>
#include <cctype>
#include <filesystem>

// Contents of this file were written with the help of ChatGPT

std::string getDataFilename(const std::string& headerFilename)
{
    std::filesystem::path path(headerFilename);
    path.replace_filename(path.stem().string() + "Data.dat");
    return path.string();
}

std::string removePath(const std::string& filename)
{
    std::filesystem::path path(filename);
    return path.filename().string();
}

std::string addPath(const std::string& srcFilename, const std::string& dstFilename)
{
    std::filesystem::path srcPath(srcFilename);
    std::filesystem::path newPath = srcPath.parent_path() / dstFilename;
    return newPath.string();
}

std::string addStepCounter(const std::string& filename, const int step)
{
    std::filesystem::path path(filename);
    std::string stem = path.stem().string();
    std::string extension = path.extension().string();
    path.replace_filename(stem + "_" + std::to_string(step) + extension);
    return path.string();
}

std::string removeStepCounter(const std::string& filename)
{
    std::filesystem::path path(filename);
    std::string stem = path.stem().string();
    std::string extension = path.extension().string();
    size_t pos = stem.find_last_of('_');
    if(pos != std::string::npos)
    {
        std::string suffix = stem.substr(pos + 1);
        if(!suffix.empty() &&
           std::all_of(suffix.begin(), suffix.end(),
                       [](unsigned char c) {return std::isdigit(c);}))
        {
            stem = stem.substr(0, pos);
        }
    }
    path.replace_filename(stem + extension);
    return path.string();
}

std::string setExtension(const std::string& filename, const std::string& extension)
{
    std::filesystem::path path(filename);
    std::string ext = extension;
    if(!ext.empty() && ext[0] != '.')
    {
        ext = "." + ext;
    }
    path.replace_extension(ext);
    return path.string();
}

void createDirForFile(const std::string& filename)
{
    std::filesystem::path path(filename);
    std::filesystem::path dir = path.parent_path();
    if(!std::filesystem::exists(dir))
    {
        std::filesystem::create_directory(dir);
    }
}
