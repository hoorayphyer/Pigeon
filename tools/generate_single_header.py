#!/usr/bin/env python3
# This file is a part of toml++ and is subject to the the terms of the MIT license.
# Copyright (c) Mark Gillard <mark.gillard@outlook.com.au>
# See https://github.com/marzer/tomlplusplus/blob/master/LICENSE for the full license text.
# SPDX-License-Identifier: MIT

import sys
import utils
import re
from pathlib import Path
from io import StringIO




class Preprocessor:

	__re_includes = re.compile(r'^\s*#\s*include\s+"(.+?)"', re.I | re.M)

	def __init__(self, file):
		self.__processed_includes = []
		self.__current_level = 0
		self.__directory_stack = [ Path.cwd() ]
		self.__entry_root = ''
		self.__string = self.__preprocess(file)

	def __preprocess(self, incl):

		if not isinstance(incl, (Path, str)): # a regex match object
			incl = incl.group(1).strip()
		if isinstance(incl, str):
			incl = Path(incl.strip().replace('\\',r'/'))
		if not incl.is_absolute():
			incl = Path(self.__directory_stack[-1], incl).resolve()
		if incl in self.__processed_includes:
			return ''
		if self.__current_level == 0 and self.__entry_root == '':
			self.__entry_root = str(incl.parent).replace('\\',r'/')

		self.__processed_includes.append(incl)
		self.__directory_stack.append(incl.parent)

		text = utils.read_all_text_from_file(incl, logger=True).strip() + '\n'
		text = text.replace('\r\n', '\n') # convert windows newlines
		self.__current_level += 1
		text = self.__re_includes.sub(lambda m : self.__preprocess(m), text, 0)
		self.__current_level -= 1

		if self.__current_level == 1:
			header = str(incl).replace('\\',r'/')
			if header.startswith(self.__entry_root):
				header = header[len(self.__entry_root):].strip('/')
			header = utils.make_divider(header, 10, pattern = r'*')
			text = f'{header}\n\n{text}'

		self.__directory_stack.pop()
		return '\n\n' + text + '\n\n'

	def __str__(self):
		return self.__string



def main():

	# establish local directories
	root_dir = utils.entry_script_dir().parent
	include_dir = Path(root_dir, 'include', 'toml++')

	# preprocess header(s)
	toml_h = str(Preprocessor(Path(include_dir, 'toml.h')))

	# strip various things:
	if 1:
		# 'pragma once'
		toml_h = re.sub(r'^\s*#\s*pragma\s+once\s*$', '', toml_h, 0, re.I | re.M)
		# trailing whitespace
		toml_h = re.sub('([^ \t])[ \t]+\n', r'\1\n', toml_h)
		# explicit 'strip this' blocks
		toml_h = re.sub(r'(?:\n[ \t]*)?//[#!][ \t]*[{][{].*?//[#!][ \t]*[}][}].*?\n', '\n', toml_h, flags=re.S)
		# spdx license identifiers
		toml_h = re.sub(r'^\s*//\s*SPDX-License-Identifier:.+?$', '', toml_h, 0, re.I | re.M)
		# magic comments
		blank_line = r'(?:[ \t]*\n)'
		comment_line = r'(?:[ \t]*//(?:[/#!<]| ?(?:---|===|\^\^\^|vvv))[^\n]*\n)'
		toml_h = re.sub(rf'\n{comment_line}{blank_line}+{comment_line}', '\n', toml_h)
		toml_h = re.sub(rf'([{{,])\s*\n(?:{comment_line}|{blank_line})+', r'\1\n', toml_h)
		toml_h = re.sub(rf'{comment_line}+', '\n', toml_h)
		# trailing whitespace
		toml_h = re.sub('([^ \t])[ \t]+\n', r'\1\n', toml_h)
		# double blank lines
		toml_h = re.sub('\n(?:[ \t]*\n[ \t]*)+\n', '\n\n', toml_h)
		# weird spacing edge case between } and pp directives
		toml_h = re.sub('\n[}]\n#', r'\n}\n\n#', toml_h, re.S)
		# blank lines following opening brackets or a comma
		toml_h = re.sub(r'([^@][({,])\n\n', r'\1\n', toml_h)
		# blank lines preceeding closing brackets
		toml_h = re.sub(r'\n\n([ \t]*[})])', r'\n\1', toml_h)
		toml_h = toml_h.strip() + '\n'

	# change TOML_LIB_SINGLE_HEADER to 1
	toml_h = re.sub(
		'#\s*define\s+TOML_LIB_SINGLE_HEADER\s+[0-9]+',
		'#define	TOML_LIB_SINGLE_HEADER 1',
		toml_h, 0, re.I
	)

	# read version number
	version_h = utils.read_all_text_from_file(Path(include_dir, 'impl/version.h'), logger=True)
	match = re.search(
			r'#\s*define\s+TOML_LIB_MAJOR\s+([0-9]+)[^0-9].*'
			+ r'#\s*define\s+TOML_LIB_MINOR\s+([0-9]+)[^0-9].*'
			+ r'#\s*define\s+TOML_LIB_PATCH\s+([0-9]+)[^0-9]',
			version_h, re.I | re.S)
	if match is None:
		raise Exception("could not find TOML_LIB_MAJOR, TOML_LIB_MINOR or TOML_LIB_PATCH impl/version.h")
	version = rf'{int(match[1])}.{int(match[2])}.{int(match[3])}'
	print(rf'Library version: {version}')

	# build the preamble (license etc)
	preamble = []
	preamble.append(rf'''
// toml++ v{version}
// https://github.com/marzer/tomlplusplus
// SPDX-License-Identifier: MIT''')
	preamble.append(r'''
// -         THIS FILE WAS ASSEMBLED FROM MULTIPLE HEADER FILES BY A SCRIPT - PLEASE DON'T EDIT IT DIRECTLY            -
//
// If you wish to submit a contribution to toml++, hooray and thanks! Before you crack on, please be aware that this
// file was assembled from a number of smaller files by a python script, and code contributions should not be made
// against it directly. You should instead make your changes in the relevant source file(s). The file names of the files
// that contributed to this header can be found at the beginnings and ends of the corresponding sections of this file.''')
	preamble.append(r'''
// TOML Language Specifications:
// latest:      https://github.com/toml-lang/toml/blob/master/README.md
// v1.0.0:      https://toml.io/en/v1.0.0
// v0.5.0:      https://toml.io/en/v0.5.0
// changelog:   https://github.com/toml-lang/toml/blob/master/CHANGELOG.md''')
	preamble.append(utils.read_all_text_from_file(Path(utils.entry_script_dir(), '..', 'LICENSE').resolve(), logger=True))

	# write the output
	with StringIO(newline='\n') as output:

		# build in a string buffer
		write = lambda txt, end='\n': print(txt, file=output, end=end)
		if (len(preamble) > 0):
			write(utils.make_divider())
		for pre in preamble:
			write('//')
			for line in pre.strip().splitlines():
				if len(line) == 0:
					write('//')
					continue
				if not line.startswith('//'):
					write('// ', end = '')
				write(line)
			write('//')
			write(utils.make_divider())
		write(toml_h)
		write('')

		output_str = output.getvalue().strip()

		# analyze the output to find any potentially missing #undefs
		if 1:
			re_define = re.compile(r'^\s*#\s*define\s+([a-zA-Z0-9_]+)(?:$|\s|\()')
			re_undef = re.compile(r'^\s*#\s*undef\s+([a-zA-Z0-9_]+)(?:$|\s|//)')
			defines = dict()
			for output_line in output_str.splitlines():
				defined = True
				m = re_define.match(output_line)
				if not m:
					defined = False
					m = re_undef.match(output_line)
				if m:
					defines[m.group(1)] = defined
			ignore_list = ( # macros that are meant to stay public (user configs etc)
				r'INCLUDE_TOMLPLUSPLUS_H',
				r'TOML_ALL_INLINE',
				r'TOML_API',
				r'TOML_CONFIG_HEADER',
				r'TOML_EXCEPTIONS',
				r'TOML_HEADER_ONLY',
				r'TOML_LANG_MAJOR',
				r'TOML_LANG_MINOR',
				r'TOML_LANG_PATCH',
				r'TOML_LARGE_FILES',
				r'TOML_LIB_MAJOR',
				r'TOML_LIB_MINOR',
				r'TOML_LIB_PATCH',
				r'TOML_LIB_SINGLE_HEADER',
				r'TOML_MAX_NESTED_VALUES',
				r'TOML_OPTIONAL_TYPE',
				r'TOML_PARSER',
				r'TOML_SMALL_FLOAT_TYPE',
				r'TOML_SMALL_INT_TYPE',
				r'TOML_UNDEF_MACROS',
				r'TOML_UNRELEASED_FEATURES',
				r'TOML_WINDOWS_COMPAT',
				r'POXY_IMPLEMENTATION_DETAIL',
			)
			set_defines = []
			for define, currently_set in defines.items():
				if currently_set and define not in ignore_list:
					set_defines.append(define)
			if len(set_defines) > 0:
				set_defines.sort()
				print(f"Potentially missing #undefs:")
				for define in set_defines:
					print(f"\t#undef {define}")

		# write the output file
		output_file_path = Path(utils.entry_script_dir(), '..', 'toml.hpp').resolve()
		print("Writing to {}".format(output_file_path))
		with open(output_file_path,'w', encoding='utf-8', newline='\n') as output_file:
			print(output_str, file=output_file)


if __name__ == '__main__':
	utils.run(main, verbose=True)
