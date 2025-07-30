import sys
import vinfo

version = vinfo.version
release = vinfo.release
patch = vinfo.patch

def get_vstr(version, release, patch):
    vstr = str(version) + '.'
    vstr += str(release) + '.'
    vstr += str(patch)
    return vstr

# run using: python vcon mode
# where mode = patch, release or version depending on the version update you
# require.

if __name__ == "__main__":

    print("")
    print("Running Version Controller")
    print("--------------------------")

    vstr_before = get_vstr(version, release, patch)

    mode = str(sys.argv[1])

    if mode == 'patch':
        patch += 1
    elif mode == 'release':
        patch = 0
        release += 1
    elif mode == 'version':
        patch = 0
        release = 0
        version += 1
    elif mode == 'test':
        pass

    print("")
    print("Mode = ", mode)

    vstr_after = get_vstr(version, release, patch)

    print("")
    print("Current version:", vstr_before)
    print("Updated version:", vstr_after)

    with open('vinfo.py', 'w') as file:
        file.write('version='+str(version)+'\n')
        file.write('release='+str(release)+'\n')
        file.write('patch='+str(patch)+'\n')
        file.write('vstr="'+vstr_after+'"\n')

    print("")
    print("Update README.md")

    readme_file = open("README.md", "r")
    list_of_lines = readme_file.readlines()
    for i in range(0, len(list_of_lines)):
        if list_of_lines[i][:17] == '| Version       |':
            list_of_lines[i] = '| Version       | ' + vstr_after
            for j in range(0, 38-len(vstr_after)):
                list_of_lines[i] += ' '
            list_of_lines[i] += '|\n'
    readme_file.close()

    readme_file = open("README.md", "w")
    readme_file.writelines(list_of_lines)
    readme_file.close()

    print("")
    print("Update index.rst")

    index_rst_file = open("docs/source/index.rst", "r")
    list_of_lines = index_rst_file.readlines()
    for i in range(0, len(list_of_lines)):
        if list_of_lines[i][:17] == '| Version       |':
            list_of_lines[i] = '| Version       | ' + vstr_after
            for j in range(0, 40-len(vstr_after)):
                list_of_lines[i] += ' '
            list_of_lines[i] += '|\n'
    index_rst_file.close()

    index_rst_file = open("docs/source/index.rst", "w")
    index_rst_file.writelines(list_of_lines)
    index_rst_file.close()

    print("")
    print("Update VERSION")

    setup_file = open("VERSION", "r")
    list_of_lines = setup_file.readlines()
    for i in range(0, 1):
        list_of_lines[i] = '%s' % vstr_after
    setup_file.close()

    setup_file = open("VERSION", "w")
    setup_file.writelines(list_of_lines)
    setup_file.close()

    print("")
