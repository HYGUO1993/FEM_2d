# -*- coding: utf-8 -*-
"""GitHub Release 附件上传工具 (直接二进制 body, 避免 multipart 污染)

用法:
  python upload_release.py <tag> <本地安装包路径> [--name 附件名]

注意:
  - 必须用直接二进制 body (Content-Type: application/octet-stream),
    multipart 上传会把文件污染成 multipart 报文 (历史教训)
  - 需本机 git credential 已登录 GitHub
"""
import subprocess, json, sys, urllib.request, urllib.error, os

REPO = "HYGUO1993/femlab-studio"


def get_token():
    r = subprocess.run(['git', 'credential', 'fill'],
                       input='protocol=https\nhost=github.com\n\n',
                       capture_output=True, text=True)
    for line in r.stdout.splitlines():
        if line.startswith('password='):
            return line.split('=', 1)[1].strip()
    return None


def api(method, url, token, data=None, headers=None):
    h = {'Authorization': 'Bearer ' + token, 'User-Agent': 'femlab-release',
         'Accept': 'application/vnd.github+json'}
    if headers:
        h.update(headers)
    req = urllib.request.Request(url, data=data, method=method, headers=h)
    return urllib.request.urlopen(req)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    tag, local_path = sys.argv[1], sys.argv[2]
    name = None
    if '--name' in sys.argv:
        name = sys.argv[sys.argv.index('--name') + 1]
    if not name:
        name = os.path.basename(local_path).replace(' ', '.')

    token = get_token()
    if not token:
        print('无法获取 GitHub token'); sys.exit(1)

    # 定位 release
    rel = json.loads(api('GET', f'https://api.github.com/repos/{REPO}/releases/tags/{tag}',
                         token).read())
    rel_id = rel['id']
    print(f'Release: {rel["name"]} (id={rel_id})')

    # 删除同名旧附件
    for a in rel['assets']:
        if a['name'] == name:
            api('DELETE', f'https://api.github.com/repos/{REPO}/releases/assets/{a["id"]}', token)
            print(f'已删除旧附件: {name}')

    # 上传 (直接二进制 body)
    data = open(local_path, 'rb').read()
    url = (f'https://uploads.github.com/repos/{REPO}/releases/{rel_id}/'
           f'assets?name={name}')
    try:
        resp = json.loads(api('POST', url, token, data=data,
                              headers={'Content-Type': 'application/octet-stream'}).read())
        print(f'上传成功: {resp["name"]} ({resp["size"]} bytes)')
        print(f'下载: {resp["browser_download_url"]}')
    except urllib.error.HTTPError as e:
        print(f'上传失败 HTTP {e.code}: {e.read().decode("utf-8", "ignore")[:500]}')
        sys.exit(1)


if __name__ == '__main__':
    main()
