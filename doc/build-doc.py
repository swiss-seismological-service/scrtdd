#!/usr/bin/env python
# -*- coding: utf-8 -*

############################################################################
# Copyright (C) gempa GmbH                                                 #
# All rights reserved.                                                     #
# Contact: gempa GmbH (seiscomp-dev@gempa.de)                              #
#                                                                          #
# GNU Affero General Public License Usage                                  #
# This file may be used under the terms of the GNU Affero                  #
# Public License version 3.0 as published by the Free Software Foundation  #
# and appearing in the file LICENSE included in the packaging of this      #
# file. Please review the following information to ensure the GNU Affero   #
# Public License version 3.0 requirements will be met:                     #
# https://www.gnu.org/licenses/agpl-3.0.html.                              #
#                                                                          #
# Other Usage                                                              #
# Alternatively, this file may be used in accordance with the terms and    #
# conditions contained in a signed written agreement between you and       #
# gempa GmbH.                                                              #
############################################################################


from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals)

import sys
import os
import glob
import shutil
import codecs
import re
import getopt
try:
    # Python 2.5
    from xml.etree import ElementTree

except ImportError:
    from elementtree import ElementTree


def copytree(src, dst, symlinks=False, ignore=None):
    names = os.listdir(src)
    if ignore is not None:
        ignored_names = ignore(src, names)
    else:
        ignored_names = set()

    try:
        os.makedirs(dst)
    except BaseException:
        pass

    errors = []
    for name in names:
        if name in ignored_names:
            continue
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                copytree(srcname, dstname, symlinks, ignore)
            else:
                # Will raise a SpecialFileError for unsupported file types
                shutil.copy2(srcname, dstname)
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except EnvironmentError as why:
            errors.append((srcname, dstname, str(why)))
    try:
        shutil.copystat(src, dst)
    except OSError as why:
        if WindowsError is not None and isinstance(why, WindowsError):
            # Copying file access times may fail on Windows
            pass
        else:
            errors.append((src, dst, str(why)))
    if errors:
        raise EnvironmentError(errors)


cmdline_templates = {}


def find_doc_dirs(directory):
    # pylint: disable=W0621
    visited = set()
    # The followlinks option has been added with Python 2.6
    if sys.version_info >= (2, 6):
        for root, _, _ in os.walk(directory, followlinks=True):
            if os.path.basename(root) == "descriptions":
                abs_root = os.path.abspath(os.path.realpath(root))
                if not abs_root in visited:
                    visited.add(abs_root)
                    yield abs_root
    else:
        for root, _, _ in os.walk(directory):
            if os.path.basename(root) == "descriptions":
                abs_root = os.path.abspath(os.path.realpath(root))
                if not abs_root in visited:
                    visited.add(abs_root)
                    yield abs_root


def escapeString(string):
    return re.sub(r"([=*\(\)|\-!@~\"&/\\\^\$\=])", r"\\\1", string)


def tagname(element):
    names = element.tag.split("}")
    if not names:
        return ""
    return names.pop()


def escape(txt):
    # return txt.replace("*", "\\*")
    return escapeString(txt)


def xml_desc_lines(n):
    desc_node = n.find('description')
    if desc_node is not None and desc_node.text is not None:
        return [l.strip()
                for l in desc_node.text.strip().replace("\r", "").split('\n')]
    return []


def xml_collect_params(param_nodes, struct_nodes, group_nodes, prefix):
    if param_nodes is None and group_nodes is None:
        return ""
    options = ""

    for param_node in param_nodes:
        name = param_node.get('name')
        type = param_node.get('type')
        unit = param_node.get('unit')

        if name is None:
            continue
        options += "\n.. confval:: %s%s\n\n" % (prefix, name)
        default = param_node.get('default')
        if type:
            options += "   Type: *%s*\n\n" % type
        if unit:
            options += "   Unit: *%s*\n\n" % unit

        desc = xml_desc_lines(param_node)

        # Description available
        if desc:
            for line in desc:
                options += "   %s\n" % escape(line)
            if default:
                options += "   Default is ``%s``." % default
        # No description, but default
        elif default:
            options += "   Default is ``%s``." % default
        # Nothing
        else:
            options += "   *No description available*"
        options += "\n"

    for struct_node in struct_nodes:
        struct_prefix = prefix + r'$name.'
        options += "\n"
        options += ".. note::\n\n"
        options += "   **%s\***\n" % struct_prefix

        desc = xml_desc_lines(struct_node)
        if desc:
            for l in desc:
                options += "   *" + l + "*\n"

        options += "   $name is a placeholder for the name to be used"
        link = struct_node.get('link')
        if link:
            options += " and needs to be added to :confval:`%s` to become active.\n\n" % link
            options += "   .. code-block:: sh\n\n"
            options += "      %s = a,b\n" % link
            options += "      %sa.value1 = ...\n" % prefix
            options += "      %sb.value1 = ...\n" % prefix
            options += "      # c is not active because it has not been added\n"
            options += "      # to the list of %s\n" % link
            options += "      %sc.value1 = ...\n" % prefix
        else:
            options += ".\n"
        options += "\n"

        options += xml_collect_params(
            struct_node.findall('parameter'),
            struct_node.findall('struct'),
            struct_node.findall('group'),
            struct_prefix)

    for group_node in group_nodes:
        group_prefix = prefix + group_node.get('name') + '.'
        desc = xml_desc_lines(group_node)
        if desc:
            options += "\n"
            options += ".. note::\n"
            options += "   **%s\***\n" % group_prefix
            for l in desc:
                options += "   *" + l + "*\n"
            options += "\n\n"

        options += xml_collect_params(
            group_node.findall('parameter'),
            group_node.findall('struct'),
            group_node.findall('group'),
            group_prefix)

    return options


def xml_collect_options(mod_node, prefix=""):
    cfg_node = mod_node.find('configuration')
    if cfg_node is None:
        return ""

    param_nodes = cfg_node.findall('parameter')
    struct_nodes = cfg_node.findall('struct')
    group_nodes = cfg_node.findall('group')

    options = xml_collect_params(
        param_nodes,
        struct_nodes,
        group_nodes,
        prefix)
    return options


def xml_collect_cmdline(mod_node, publicids_only):
    cmd_node = mod_node.find("command-line")
    if cmd_node is None:
        return ""

    group_nodes = cmd_node.findall("group")
    if not group_nodes:
        return ""

    options = '''

Command Line
============

'''

    synopsis = cmd_node.find("synopsis")
    if synopsis is not None:
        print("Found synopsis")
        options += ":program:`%s`\n\n" % synopsis.text.strip()

    desc = xml_desc_lines(cmd_node)
    if desc:
        for line in desc:
            options += "%s\n" % escape(line)
        options += "\n"

    for group_node in group_nodes:
        name = group_node.get('name')
        if name:
            options += '''
%s
%s

''' % (name, '-' * len(name))

        if not publicids_only:
            optionref_nodes = group_node.findall('optionReference')
            for optionref_node in optionref_nodes:
                try:
                    publicID = optionref_node.text.strip()
                except BaseException:
                    sys.stderr.write("WARNING: publicID is empty\n")
                    continue
                if publicID not in cmdline_templates:
                    sys.stderr.write(
                        "WARNING: option with publicID '%s' is not available\n" %
                        publicID)
                    continue
                options += cmdline_templates[publicID]

        option_nodes = group_node.findall('option')
        for option_node in option_nodes:
            flag = option_node.get('flag')
            long_flag = option_node.get('long-flag')
            param_ref = option_node.get('param-ref')
            flags = ""
            if flag:
                flags += "-" + flag
            if long_flag:
                if flags:
                    flags += ", "
                flags += "--" + long_flag

            arg = option_node.get('argument')
            if arg:
                arg = " " + arg
            else:
                arg = ""

            option = ".. option:: %s%s\n\n" % (flags, arg)

            if param_ref:
                option += "   Overrides configuration parameter :confval:`%s`.\n" % param_ref
            desc = xml_desc_lines(option_node)
            if desc:
                for line in desc:
                    option += "   %s\n" % escape(line)
            option += "\n"

            publicID = option_node.get('publicID')
            if publicID:
                cmdline_templates[publicID] = option

            if not publicids_only:
                options += option

    return options

# ------------------------------------------------------------------------------


def resolveVariables(inFile, outFile, d):
    print("Generating {}".format(outFile))
    with open(inFile, 'r') as fIn:
        t = fIn.read()
        with open(outFile, 'w') as fOut:
            for k, v in d.items():
                search = '${generator.' + k + '}'
                # print " replacing {} with {}".format(search, v)
                t = t.replace(search, v)

            fOut.write(t)


# ------------------------------------------------------------------------------
def usage():
    print('''
usage: %s [OPTIONS] sourceDir [outputDir]

Options:
  -h,--help      print this help text
  -v,--version   version string to use in replace in conf.py template
  --sc3          build for SeisComP3
  --skip-category do not generate module category menu
''' % sys.argv[0])


# ------------------------------------------------------------------------------
try:
    opts, args = getopt.getopt(sys.argv[1:], "hv:", [
                               "help", "version=", "sc3", "skip-category"])
except getopt.GetoptError as err:
    # print help information and exit:
    print(str(err))  # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

# key value pairs to resolve in templates
target_sc3 = False
skip_category = False
resolveDict = {}
for o, a in opts:
    if o in ('-h', '--help'):
        usage()
        sys.exit(0)
    elif o in ('-v', '--version'):
        resolveDict['param.version'] = a
    elif o in ('--sc3'):
        target_sc3 = True
    elif o in ('--skip-category'):
        skip_category = True
    else:
        sys.stderr.write('Unknown option: {}\n'.format(o))
        usage()
        sys.exit(1)

if len(args) < 1:
    sys.stderr.write("No source directory given\n")
    usage()
    sys.exit(1)

exe_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
src_dir = args[0]
sys.stderr.write("Building configuration from '%s'\n" % src_dir)

if len(args) > 1:
    out_build_dir = args[1]
else:
    out_build_dir = "build-doc"


toc_tmp_dir = "toc"

out_dir = os.path.join(out_build_dir, "src")
out_struct_dir = os.path.join(out_dir, toc_tmp_dir)

if os.path.exists(out_dir):
    print("Cleaning target dir")
    shutil.rmtree(out_dir)


if not os.path.exists(out_struct_dir):
    try:
        os.makedirs(out_struct_dir)
    except OSError as e:
        sys.stderr.write(
            "ERROR: creating build directory '%s' failed: %s\n" %
            (out_struct_dir, str(e)))
        sys.exit(1)

apps_doc_dir = os.path.join(src_dir, "apps")
apps_doc_base_dir = os.path.join(exe_dir, "apps")

categories = {}
app_nodes = []
app_plugin_nodes = {}
app_binding_nodes = {}
global_node = None

doc_dirs = list(find_doc_dirs(".."))
doc_dirs.append(apps_doc_dir)
doc_dirs.append(apps_doc_base_dir)

# Read all apps/*.xml files and grab the "module" node
# Plugins and bindings will follow latera
for doc_dir in doc_dirs:
    for f in glob.glob(doc_dir + "/*.xml"):
        try:
            xmlTree = ElementTree.parse(f)
        except ElementTree.ParseError as e:
            sys.stderr.write("ERROR: %s: parsing failed: %s\n" % (f, str(e)))
            sys.exit(1)
        except IOError as e:
            sys.stderr.write("Warning: %s\n" % str(e))
            continue
    
        root = xmlTree.getroot()
        if root is None:
            sys.stderr.write("ERROR: %s: no root tag defined\n" % f)
            sys.exit(1)
    
        if tagname(root) != "seiscomp":
            sys.stderr.write(
                "ERROR: %s: invalid root tag: expected 'seiscomp'\n" %
                f)
            sys.exit(1)
    
        for node in root:
            tag = tagname(node)
            if tag == "module":
                name = node.get('name')
                if not name:
                    sys.stderr.write("ERROR: %s: module name not defined\n" % f)
                    sys.exit(1)
    
                # The global configuration will be handled differently
                if name == "global":
                    global_node = node
                    continue
    
                category = node.get('category')
                if not category:
                    category = "Modules"
    
                toks = category.split('/')
                ct = []
                for t in toks[:-1]:
                    ct = ct + [t]
                    if '/'.join(ct) not in categories:
                        categories['/'.join(ct)] = []
    
                # Add node to category map
                try:
                    categories[category].append(node)
                except BaseException:
                    categories[category] = [node]
    
                app_nodes.append(node)
    
            elif tag == "plugin":
                name = node.get('name')
                if not name:
                    sys.stderr.write("ERROR: %s: plugin name not defined\n" % f)
                    sys.exit(1)
    
                try:
                    extends = node.find('extends').text.strip()
                except BaseException:
                    sys.stderr.write("ERROR: %s: plugin extends not defined\n" % f)
                    sys.exit(1)
    
                # Save all plugins of a module
                try:
                    app_plugin_nodes[extends].append(node)
                except BaseException:
                    app_plugin_nodes[extends] = [node]
    
            elif tag == "binding":
                modname = node.get('module')
                if not modname:
                    sys.stderr.write("ERROR: %s: binding module not defined\n" % f)
                    sys.exit(1)
    
                if node.get('category') and not node.get('name'):
                    sys.stderr.write(
                        "ERROR: %s: binding category defined but no name given\n" %
                        f)
                    sys.exit(1)
    
                # Save all plugins of a module
                try:
                    app_binding_nodes[modname].append(node)
                except BaseException:
                    app_binding_nodes[modname] = [node]


# categories hold pair [category, [node1, node2, ...]]
# app_refs holds [reference_link, [section name, nodes]]
app_refs = {}

for key, value in sorted(categories.items()):
    dirs = key.split('/')
    section = dirs[-1]
    section_path = dirs[:-1]

    if section_path:
        link = os.path.join(toc_tmp_dir, *section_path).lower()
        path = os.path.join(out_dir, link)
        if not os.path.exists(path):
            try:
                os.makedirs(path)
            except OSError as e:
                sys.stderr.write(
                    "ERROR: creating path '%s' failed: %s\n" %
                    (path, str(e)))
                sys.exit(1)
        link = os.path.join(link, section.lower())
    else:
        link = os.path.join(toc_tmp_dir, section.lower())

    app_refs[link] = [section_path, section, value]


placeholder_app_refs = ""

print("Generating document structure")

for ref, nodes in sorted(app_refs.items()):
    section_path = nodes[0]
    section = nodes[1]
    xml_nodes = nodes[2]

    if not section_path:
        placeholder_app_refs += "   /" + ref + "\n"

    if skip_category:
        continue

    try:
        f = open(os.path.join(out_dir, ref) + ".rst", "w")
    except BaseException:
        sys.stderr.write("ERROR: unable to create index file: %s\n" % ref)
        sys.exit(1)

    # Create section files
    f.write('#' * len(section) + "\n")
    f.write("%s\n" % section)
    f.write('#' * len(section) + "\n")
    f.write("\n")
    f.write(".. toctree::\n")
    f.write("   :maxdepth: 2\n")
    f.write("\n")

    for n in sorted(xml_nodes, key=lambda n: n.get('name')):
        f.write("   /apps/" + n.get('name').lower() + "\n")

    my_path = section_path + [section]
    child_refs = []
    for ref2, nodes2 in app_refs.items():
        if nodes2[0] == my_path:
            child_refs.append(ref2)

    for cref in sorted(child_refs):
        f.write("   /" + cref + "\n")

    f.close()

# Prepend global
if global_node is not None:
    placeholder_app_refs = "   /apps/global\n" + placeholder_app_refs

resolveDict['refs.apps'] = placeholder_app_refs

# Create conf.py
resolveVariables(os.path.join(src_dir, "templates", "conf.py"),
                 os.path.join(out_dir, "conf.py"), resolveDict)

# Copy base directory
print("Copy base directory")
out_base_dir = os.path.join(out_dir, "base")
copytree(os.path.join(exe_dir, "base"), out_base_dir)
copytree(os.path.join(src_dir, "base"), out_base_dir)

out_media_dir = os.path.join(out_dir, "apps", "media")

# Copy base media files
if os.path.exists(os.path.join(exe_dir, "apps", "media")):
    print("Copy base app media files")
    copytree(os.path.join(exe_dir, "apps", "media"), out_media_dir)

# Copy media files
if os.path.exists(os.path.join(src_dir, "apps", "media")):
    print("Copy media files")
    copytree(os.path.join(src_dir, "apps", "media"), out_media_dir)

# Create index.rst
resolveVariables(os.path.join(src_dir, "templates", "index.rst"),
                 os.path.join(out_dir, "index.rst"), resolveDict)

# Create application .rst files
print("Generating app .rst files")
app_path = os.path.join(out_dir, "apps")
if not os.path.exists(app_path):
    try:
        os.makedirs(app_path)
    except OSError as e:
        sys.stderr.write(
            "ERROR: creating path '%s' failed: %s\n" %
            (app_path, str(e)))
        sys.exit(1)


if global_node is None:
    node_list = app_nodes
else:
    node_list = app_nodes + [global_node]

# First pass, collect commandline templates for options with a publicID
for n in node_list:
    xml_collect_cmdline(n, True)

man_pages = []

root_dir = "seiscomp"
if target_sc3:
    root_dir = "seiscomp3"

for n in node_list:
    app_name = n.get('name')
    filename = os.path.join(app_path, (app_name + ".rst").lower())
    man_pages.append((app_name.lower(), n.get('author')))
    try:
        f = codecs.open(filename, "w", "utf-8")
    except UnicodeDecodeError as e:
        sys.stderr.write(
            "ERROR: unable to create app rst '%s': %s\n" % (filename, str(e)))
        sys.exit(1)

    desc = "\n".join(xml_desc_lines(n))
    if not desc:
        desc = app_name

    if os.path.exists(os.path.join(apps_doc_dir, app_name + ".rst")):
        doc = codecs.open(
            os.path.join(
                apps_doc_dir,
                app_name +
                ".rst"),
            "r",
            "utf-8").read()
        # Copy original .rst to .doc
        # shutil.copyfile(os.path.join(src_dir, "apps", app_name + ".rst"),
        # os.path.join(out_dir, "apps", app_name + ".doc"))
    elif os.path.exists(os.path.join(apps_doc_base_dir, app_name + ".rst")):
        doc = codecs.open(
            os.path.join(
                apps_doc_base_dir,
                app_name + ".rst"),
            "r",
            "utf-8").read()
    else:
        doc = None

    standalone = n.get('standalone')

    if app_name != "global":
        f.write('''\
.. highlight:: rst

.. _%s:

%s
%s
%s

**%s**\n''' % (app_name, "#" * len(app_name), app_name, "#" * len(app_name), desc))

        if doc is not None:
            f.write('''

Description
===========

%s
''' % doc)
    else:
        if doc:
            f.write(".. _%s:\n\n" % app_name)
            f.write(doc)

    options = xml_collect_options(n)
    plugins = app_plugin_nodes.get(app_name)
    if plugins:
        for p in plugins:
            name = p.get('name')
            desc_node = p.find('description')
            if desc_node is not None and desc_node.text is not None:
                desc = [
                    l.strip() for l in desc_node.text.strip().replace(
                        "\r", "").split('\n')]
            else:
                desc = []
            poptions = xml_collect_options(p)
            if poptions:
                options += "\n.. _%s/%s:\n\n" % (app_name, name)
                title = name + " plugin"
                options += '''
%s
%s

''' % (title, "-" * len(title))
                options += "\n".join(desc) + "\n\n"
                options += poptions

    if standalone and standalone.lower() == "true":
        if options:
            note = '''

.. note::

   %s is a standalone module and does not inherit `global options <https://docs.gempa.de/seiscomp/current/apps/global.html>`_.

''' % app_name
            cfgs = '''\
| :file:`etc/defaults/%s.cfg`
| :file:`etc/%s.cfg`
| :file:`~/.%s/%s.cfg`

''' % (app_name, app_name, root_dir, app_name)
        else:
            note = ""
            cfgs = ""
    else:
        note = ""
        if app_name != "global":
            cfgs = '''\
| :file:`etc/defaults/global.cfg`
| :file:`etc/defaults/%s.cfg`
| :file:`etc/global.cfg`
| :file:`etc/%s.cfg`
| :file:`~/.%s/global.cfg`
| :file:`~/.%s/%s.cfg`

%s inherits `global options <https://docs.gempa.de/seiscomp/current/apps/global.html>`_.

''' % (app_name, app_name, root_dir, root_dir, app_name, app_name)
        else:
            cfgs = ""

    if options or note or cfgs:
        f.write('''

Configuration
=============
%s
%s
''' % (note, cfgs))

    if options:
        f.write(options)

        bindings = app_binding_nodes.get(app_name)
        bindings_options = ""
        if bindings:
            category_map = {}
            for b in bindings:
                category = b.get('category')
                if category:
                    # Just remember category bindings and add them afterwards
                    try:
                        category_map[category].append(b)
                    except BaseException:
                        category_map[category] = [b]
                    continue
                bindings_options += xml_collect_options(b)

            for cat, bindings in sorted(category_map.items()):
                names = []
                bindings = sorted(bindings, key=lambda n: n.get('name'))
                for b in bindings:
                    names.append(b.get('name'))

                if not names:
                    continue

                bindings_options += ".. confval:: %s\n\n" % cat
                bindings_options += "   Type: *list:string*\n\n"
                bindings_options += "   Defines a list of plugin bindings to be used.\n"
                bindings_options += "   Each binding can then be configured individually.\n\n"
                bindings_options += "   Available identifiers: %s\n\n" % " ".join(
                    [":ref:`%s-%s-%s-label`" % (app_name, cat, name) for name in names])
                bindings_options += "   .. code-block:: sh\n\n"
                if len(names) >= 2:
                    bindings_options += "      # param1 and param2 are just placeholders.\n"
                    bindings_options += "      %s = %s, %s\n" % (
                        cat, names[0], names[1])
                    bindings_options += "      %s.%s1.param1 = value11\n" % (
                        cat, names[0])
                    bindings_options += "      %s.%s1.param2 = value12\n" % (
                        cat, names[0])
                    bindings_options += "      %s.%s2.param1 = value21\n" % (
                        cat, names[1])
                    bindings_options += "      %s.%s2.param2 = value22\n" % (
                        cat, names[1])
                    bindings_options += "\n"
                bindings_options += "      # To use the same binding twice, \
                                           # aliases must be used.\n"
                bindings_options += "      # Aliases are created by prepending\n \
                                           # a unique name followed by a colon\n"
                bindings_options += "      %s = %s, %s_2:%s\n" % (
                    cat, names[0], names[0], names[0])
                bindings_options += "      %s.%s.param1 = value11\n" % (
                    cat, names[0])
                bindings_options += "      %s.%s.param2 = value12\n" % (
                    cat, names[0])
                bindings_options += "      %s.%s_2.param1 = value21\n" % (
                    cat, names[0])
                bindings_options += "      %s.%s_2.param2 = value22\n" % (
                    cat, names[0])
                bindings_options += "\n"

                for b in bindings:
                    name = b.get('name')
                    bindings_options += "\n.. _%s-%s-%s-label:\n\n" % (
                        app_name, cat, name)
                    bindings_options += "%s\n" % name
                    bindings_options += "%s\n\n" % ('-' * len(name))
                    desc = xml_desc_lines(b)
                    if desc:
                        bindings_options += "\n".join(desc)
                        bindings_options += "\n"
                    bindings_options += xml_collect_options(
                        b, cat + "." + name + ".")

        if bindings_options:
            f.write('''

Bindings
========

%s
''' % bindings_options)

    f.write(xml_collect_cmdline(n, False))


# Generate plugins.doc to give an overview over available plugin's
if app_plugin_nodes:
    filename = os.path.join(out_base_dir, "plugins.doc")
    try:
        f = codecs.open(filename, "w", "utf-8")
    except UnicodeDecodeError as e:
        sys.stderr.write("ERROR: unable to create _plugins.rst: %s\n" % str(e))
        sys.exit(1)

    table = []

    for key, nodes in app_plugin_nodes.items():
        mod = ":ref:`%s<%s>`" % (key, key)
        plugs = []
        for node in nodes:
            plugs.append(":ref:`%s<%s/%s>`" %
                         (node.get('name'), key, node.get('name')))
        table.append([mod, " ".join(plugs)])

    col1_width = len("Module")
    col2_width = len("Plugin's")
    for row in table:
        if len(row[0]) > col1_width:
            col1_width = len(row[0])
        if len(row[1]) > col2_width:
            col2_width = len(row[1])

    f.write("%s  %s\n" % ("=" * col1_width, "=" * col2_width))
    f.write("%*s  %*s\n" % (col1_width, "Module", col2_width, "Plugin's"))
    f.write("%s  %s\n" % ("=" * col1_width, "=" * col2_width))
    for row in table:
        f.write("%*s  %*s\n" % (col1_width, row[0], col2_width, row[1]))
    f.write("%s  %s\n" % ("=" * col1_width, "=" * col2_width))
f = open(os.path.join(out_dir, "conf.py"), "a")
print("# -- Options for manual page output --------------------------------------------", file=f)
print("", file=f)
print("# One entry per manual page. List of tuples", file=f)
print("# (source start file, name, description, authors, manual section).", file=f)
print("man_pages = [", file=f)
for man_page in man_pages:
    author = "GFZ Potsdam"
    if man_page[1]:
        author = man_page[1]
    print("    ('apps/%s', '%s', project + u' Documentation',\
          [u'%s'], 1)," % (man_page[0], man_page[0], author), file=f)
print("]", file=f)
print("", file=f)
f.close()

print("Clean-up HTML build dir")
try:
    shutil.rmtree(os.path.join(out_build_dir, "html"))
except BaseException:
    pass

print("Clean-up MAN build dir")
try:
    shutil.rmtree(os.path.join(out_build_dir, "man1"))
except BaseException:
    pass

try:
    os.environ["PYTHONPATH"] = os.path.abspath(
        '../src/system/libs/python') + ":" + os.environ["PYTHONPATH"]
except BaseException:
    os.environ["PYTHONPATH"] = os.path.abspath('../src/system/libs/python')

os.environ["PYTHONPATH"] = os.path.join(
    exe_dir, "libs") + ":" + os.environ["PYTHONPATH"]

os.system("sphinx-build -b html %s %s" %
          (out_dir, os.path.join(out_build_dir, "html")))
#os.system("sphinx-build -b pdf %s %s" % (out_dir, os.path.join(out_build_dir, "pdf")))
os.system("sphinx-build -b man %s %s" %
          (out_dir, os.path.join(out_build_dir, "man1")))

print("Clean-up temporary files")
try:
    shutil.rmtree(os.path.join(out_build_dir, "man1", ".doctrees"))
except BaseException:
    pass

try:
    shutil.rmtree(os.path.join(out_build_dir, "html", ".doctrees"))
except BaseException:
    pass

try:
    shutil.rmtree(os.path.join(out_build_dir, "pdf", ".doctrees"))
except BaseException:
    pass
