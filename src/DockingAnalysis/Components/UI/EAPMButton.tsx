import { ButtonHTMLAttributes } from "react";

export default function EAPMButton(
  props: ButtonHTMLAttributes<HTMLButtonElement> & {
    setClass?: boolean;
    buttonColor?: string;
  }
) {
  const buttonColor = props.buttonColor ?? "blue";

  let className = `shadow-lg shadow-${buttonColor}-500/50 font-medium rounded-lg text-sm px-5 py-2.5 text-center me-2 mb-2 ${props.className}`;

  if (props.setClass === true || props.setClass === undefined) {
    className =
      `border border-${buttonColor}-500 text-black bg-gradient-to-r from-${buttonColor}-500 via-${buttonColor}-600 to-${buttonColor}-700 hover:bg-gradient-to-br ` +
      className;
  }

  return (
    <button {...props} className={className}>
      {props.children}
    </button>
  );
}
