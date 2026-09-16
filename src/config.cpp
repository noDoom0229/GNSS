/*-----------------------------------------------------------------------------
    config.cpp  —— 读取 INI 风格配置文件（对应报告 3.2.1 config.cpp）

    每行一个配置项：key = value
    以 # 或 ; 开头的行是注释，空行跳过。
-----------------------------------------------------------------------------*/
#include "SPP.h"
#include <map>
#include <cstdlib>

// 去掉字符串两端的空白字符
static string trim(const string& s)
{
    size_t begin = s.find_first_not_of(" \t\r\n");
    if (begin == string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(begin, end - begin + 1);
}

// 从键值表中取字符串 / 整数 / 浮点数，取不到时保持默认值
static string get_str(map<string, string>& kv, const string& key, const string& def)
{
    if (kv.count(key)) return kv[key];
    return def;
}

static int get_int(map<string, string>& kv, const string& key, int def)
{
    if (kv.count(key)) return atoi(kv[key].c_str());
    return def;
}

static double get_double(map<string, string>& kv, const string& key, double def)
{
    if (kv.count(key)) return atof(kv[key].c_str());
    return def;
}

bool load_config(const string& filename, ConfigInfo& config)
{
    ifstream fin(filename.c_str());
    if (!fin.is_open())
    {
        printf("Cannot open config file: %s\n", filename.c_str());
        return false;
    }

    map<string, string> kv;
    string line;
    while (getline(fin, line))
    {
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] == '#' || line[0] == ';') continue;

        size_t pos = line.find('=');
        if (pos == string::npos) continue;

        string key = trim(line.substr(0, pos));
        string value = trim(line.substr(pos + 1));
        if (key.empty()) continue;
        kv[key] = value;
    }
    fin.close();

    config.BasBinFile = get_str(kv, "BasBinFile", "");
    config.RovBinFile = get_str(kv, "RovBinFile", "");
    config.BasIP = get_str(kv, "BasIP", "");
    config.RovIP = get_str(kv, "RovIP", "");
    config.BasPort = get_int(kv, "BasPort", 0);
    config.RovPort = get_int(kv, "RovPort", 0);
    config.NMEAOutputFile = get_str(kv, "NMEAOutputFile", "output/result.nmea");
    config.OutputFile = get_str(kv, "OutputFile", "output/result.txt");
    config.DecodeOutput = get_int(kv, "DecodeOutput", 0);
    config.BasObsFile = get_str(kv, "BasObsFile", "output/base_obs.txt");
    config.RovObsFile = get_str(kv, "RovObsFile", "output/rove_obs.txt");
    config.NavOutFile = get_str(kv, "NavOutFile", "output/nav.txt");
    config.PosMode = get_int(kv, "PosMode", 0);
    config.CalcMode = get_int(kv, "CalcMode", 0);
    config.BasX = get_double(kv, "BasX", 0.0);
    config.BasY = get_double(kv, "BasY", 0.0);
    config.BasZ = get_double(kv, "BasZ", 0.0);
    config.RovX = get_double(kv, "RovX", 0.0);
    config.RovY = get_double(kv, "RovY", 0.0);
    config.RovZ = get_double(kv, "RovZ", 0.0);
    config.ElevThreshold = get_double(kv, "ElevThreshold", 10.0);
    config.PseuThreshold = get_double(kv, "PseuThreshold", 10.0);
    config.RatioThreshold = get_double(kv, "RatioThreshold", 3.0);

    // 基本检查
    if (config.PosMode == 0 && (config.BasBinFile.empty() || config.RovBinFile.empty()))
    {
        printf("Config error: PosMode = 0 requires BasBinFile and RovBinFile.\n");
        return false;
    }
    if (config.PosMode == 1 && (config.BasIP.empty() || config.RovIP.empty() ||
        config.BasPort <= 0 || config.RovPort <= 0))
    {
        printf("Config error: PosMode = 1 requires BasIP/RovIP/BasPort/RovPort.\n");
        return false;
    }
    if (config.RatioThreshold <= 0.0)
    {
        printf("Config error: RatioThreshold must be > 0.\n");
        return false;
    }
    return true;
}
